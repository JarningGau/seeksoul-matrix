#!/usr/bin/env python3
"""Generate and optionally submit seeksoul-matrix workflow commands."""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
from pathlib import Path

from _version import __version__
import workflow_input_checks as wic

STAGE_SEQUENCE = ["fastp_split", "demux_extract_bc"]
STAGE_CHOICES = [*STAGE_SEQUENCE, "all"]
SLURM_NEST_STAGE_KEYS = frozenset(STAGE_SEQUENCE)
STAGE_REQUIRED_FIELDS = {
    "fastp_split": ["r1", "r2", "number_of_split_parts"],
    "demux_extract_bc": ["barcode_whitelist"],
}
DEFAULT_BARCODE_WHITELIST = "whitelist/DD-MET5/U3CB_methylation.txt.gz"


def resolve_env_executable(name: str) -> str:
    candidate = Path(sys.executable).resolve().parent / name
    if candidate.is_file():
        return str(candidate)
    return name


def normalize_executable_setting(value: str | None, default_name: str) -> str:
    if not value or value == default_name:
        return resolve_env_executable(default_name)
    return value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate executable command scripts for seeksoul-matrix workflow."
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    parser.add_argument(
        "--workflow-config",
        help="JSON config path for workflow/sample settings.",
    )
    parser.add_argument(
        "--runner",
        choices=["local", "slurm"],
        help="Command target: local shell or slurm sbatch.",
    )
    parser.add_argument(
        "--stage",
        choices=STAGE_CHOICES,
        help="Workflow stage to generate command script for. Default: fastp_split.",
    )
    parser.add_argument("--sample-id", help="Sample identifier.")
    parser.add_argument("--r1", help="Input R1 FASTQ(.gz).")
    parser.add_argument("--r2", help="Input R2 FASTQ(.gz).")
    parser.add_argument("--work-root", help="Work root directory. Default: work.")
    parser.add_argument(
        "--fastp-threads",
        type=int,
        help="Thread count for fastp split. Default: 8.",
    )
    parser.add_argument(
        "--number-of-split-parts",
        type=int,
        help="Value passed to fastp --split.",
    )
    parser.add_argument(
        "--fastp-bin",
        help=(
            "fastp executable path or command name. "
            "Default: fastp from current Python env if available, else fastp."
        ),
    )
    parser.add_argument(
        "--barcode-whitelist",
        help=f"Cell barcode whitelist for demux. Default: {DEFAULT_BARCODE_WHITELIST}.",
    )
    parser.add_argument(
        "--barcode-hamming-distance",
        type=int,
        help="Hamming distance for demux barcode correction. Default: 1.",
    )
    parser.add_argument(
        "--gzip-level",
        type=int,
        help="gzip level for demux output FASTQ. Default: 6.",
    )
    parser.add_argument(
        "--submit",
        action="store_true",
        help="Submit immediately after generating command file.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print command and output path without writing files.",
    )
    parser.add_argument(
        "--skip-workdir-input-checks",
        action="store_true",
        help=(
            "Do not require prior-stage outputs under the sample work directory. "
            "For fastp_split, still skips r1/r2 existence checks when set. "
            "When generating --stage all, this is passed to each per-stage subprocess automatically."
        ),
    )
    parser.add_argument("--slurm-partition")
    parser.add_argument("--slurm-mem")
    parser.add_argument("--slurm-cpus-per-task", type=int)
    parser.add_argument("--slurm-output")
    parser.add_argument("--slurm-error")
    return parser.parse_args()


def quoted(parts: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in parts)


def build_fastp_split_command(args: argparse.Namespace, sample_work: Path) -> str:
    script_path = Path("scripts/fastp_split.py")
    command = [
        sys.executable,
        str(script_path),
        "--r1",
        args.r1,
        "--r2",
        args.r2,
        "--work-path",
        str(sample_work),
        "--fastp-threads",
        str(args.fastp_threads),
        "--number-of-split-parts",
        str(args.number_of_split_parts),
        "--fastp-bin",
        args.fastp_bin,
    ]
    return quoted(command)


def build_demux_chunk_command(
    args: argparse.Namespace, r1_path: Path, r2_path: Path, out_prefix: Path
) -> str:
    return quoted(
        [
            sys.executable,
            "scripts/demux_extract_bc.py",
            str(r1_path),
            str(r2_path),
            "--barcode-whitelist",
            args.barcode_whitelist,
            "--output-prefix",
            str(out_prefix),
            "--barcode-hamming-distance",
            str(args.barcode_hamming_distance),
            "--gzip-level",
            str(args.gzip_level),
        ]
    )


def build_aggregate_ct_command(demux_dir: Path) -> str:
    return quoted(
        [
            sys.executable,
            "scripts/aggregate_ct_qc.py",
            "--demux-dir",
            str(demux_dir),
        ]
    )


def build_demux_local_batch_command(
    args: argparse.Namespace, sample_work: Path, number_of_split_parts: int | None = None
) -> str:
    demux_dir = sample_work / "demux"
    aggregate_cmd = build_aggregate_ct_command(demux_dir)
    lines = [
        f"demux_dir={shlex.quote(str(demux_dir))}",
        'mkdir -p "$demux_dir"',
        "",
    ]
    chunks = discover_demux_chunks(sample_work, number_of_split_parts)
    total = len(chunks)
    for idx, (chunk_id, r1_path, r2_path, out_prefix) in enumerate(chunks, start=1):
        percent = idx * 100 // total if total else 100
        lines.append(
            f"echo '[demux] {percent:3d}% ({idx}/{total}) {chunk_id}'"
        )
        lines.append(
            build_demux_chunk_command(args, r1_path, r2_path, out_prefix)
        )
        lines.append("")
    lines.extend([aggregate_cmd, "", 'echo "[demux] done"'])
    return "\n".join(lines)


def build_demux_slurm_submit_command(sample_work: Path) -> str:
    command_dir = sample_work / "commands"
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        f'SCRIPT_DIR={shlex.quote(str(command_dir))}',
        'job_ids=""',
        'for script in "$SCRIPT_DIR"/02_demux_extract_bc_*.sbatch; do',
        '  [ -f "$script" ] || continue',
        '  jid=$(sbatch --parsable "$script")',
        '  job_ids="${job_ids}:${jid}"',
        "done",
        'if [ -z "$job_ids" ]; then',
        '  echo "[demux] no chunk sbatch scripts found"',
        "  exit 1",
        "fi",
        'sbatch --dependency=afterok"${job_ids}" "$SCRIPT_DIR/02_aggregate_ct_qc.sbatch"',
        'echo "[demux] submitted aggregate job with dependency afterok${job_ids}"',
    ]
    return "\n".join(lines)


def discover_demux_chunks(
    sample_work: Path, number_of_split_parts: int | None = None
) -> list[tuple[str, Path, Path, Path]]:
    shard_dir = sample_work / "shard_fastq"
    demux_dir = sample_work / "demux"
    shards = wic.discover_fastp_shards(shard_dir)
    if not shards and number_of_split_parts:
        shards = wic.plan_fastp_shards(shard_dir, number_of_split_parts)
    chunks: list[tuple[str, Path, Path, Path]] = []
    for chunk_id, r1_path, r2_path in shards:
        out_prefix = demux_dir / chunk_id
        chunks.append((chunk_id, r1_path, r2_path, out_prefix))
    return chunks


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def load_workflow_config(path: str) -> dict:
    config_path = Path(path)
    with config_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError("workflow config must be a JSON object")
    return data


def pick(cli_value, cfg_value):
    return cli_value if cli_value is not None else cfg_value


def validate_required_for_stage(stage: str, settings: dict) -> None:
    required = ["runner", "sample_id", *STAGE_REQUIRED_FIELDS[stage]]
    missing = [key for key in required if settings.get(key) in (None, "")]
    if missing:
        raise ValueError(f"missing required settings: {', '.join(missing)}")


def validate_inputs_for_stage(
    stage: str,
    settings: dict,
    sample_work: Path,
    *,
    skip_workdir_inputs: bool = False,
) -> None:
    if skip_workdir_inputs:
        return
    if stage == "fastp_split":
        wic.require_file("r1", wic.resolve_config_path(settings["r1"]))
        wic.require_file("r2", wic.resolve_config_path(settings["r2"]))
        wic.require_optional_executable_path("fastp_bin", settings["fastp_bin"])
    elif stage == "demux_extract_bc":
        wic.require_file(
            "barcode_whitelist",
            wic.resolve_config_path(settings["barcode_whitelist"]),
        )
        if not skip_workdir_inputs:
            shard_dir = sample_work / "shard_fastq"
            shards = wic.discover_fastp_shards(shard_dir)
            if not shards:
                raise ValueError(f"no fastp shards found under {shard_dir}")
            for _chunk_id, r1_path, r2_path in shards:
                wic.require_file(f"shard_fastq/{r1_path.name}", r1_path)
                wic.require_file(f"shard_fastq/{r2_path.name}", r2_path)
    else:
        raise ValueError(f"unsupported stage for input validation: {stage}")


def select_stage_slurm_cfg(slurm_cfg_raw: dict, stage: str) -> dict:
    if any(key in slurm_cfg_raw for key in SLURM_NEST_STAGE_KEYS):
        stage_slurm_cfg = slurm_cfg_raw.get(stage, {})
    else:
        stage_slurm_cfg = slurm_cfg_raw
    if stage_slurm_cfg is None:
        stage_slurm_cfg = {}
    if not isinstance(stage_slurm_cfg, dict):
        raise ValueError("selected slurm config must be an object")
    return stage_slurm_cfg


def resolve_settings(args: argparse.Namespace) -> dict:
    cfg: dict = {}
    if args.workflow_config:
        cfg = load_workflow_config(args.workflow_config)

    slurm_cfg_raw = cfg.get("slurm", {})
    if slurm_cfg_raw is None:
        slurm_cfg_raw = {}
    if not isinstance(slurm_cfg_raw, dict):
        raise ValueError("workflow config key 'slurm' must be an object")

    stage = pick(args.stage, cfg.get("stage")) or "fastp_split"
    if stage not in STAGE_CHOICES:
        raise ValueError(f"unsupported stage: {stage}")

    stage_slurm_cfg = select_stage_slurm_cfg(slurm_cfg_raw, stage)
    settings = {
        "runner": pick(args.runner, cfg.get("runner")),
        "stage": stage,
        "sample_id": pick(args.sample_id, cfg.get("sample_id")),
        "r1": pick(args.r1, cfg.get("r1")),
        "r2": pick(args.r2, cfg.get("r2")),
        "work_root": pick(args.work_root, cfg.get("work_root")),
        "fastp_threads": pick(args.fastp_threads, cfg.get("fastp_threads")),
        "number_of_split_parts": pick(
            args.number_of_split_parts, cfg.get("number_of_split_parts")
        ),
        "fastp_bin": pick(args.fastp_bin, cfg.get("fastp_bin")),
        "barcode_whitelist": pick(args.barcode_whitelist, cfg.get("barcode_whitelist")),
        "barcode_hamming_distance": pick(
            args.barcode_hamming_distance, cfg.get("barcode_hamming_distance")
        ),
        "gzip_level": pick(args.gzip_level, cfg.get("gzip_level")),
        "slurm_partition": pick(args.slurm_partition, stage_slurm_cfg.get("partition")),
        "slurm_mem": pick(args.slurm_mem, stage_slurm_cfg.get("mem")),
        "slurm_cpus_per_task": pick(
            args.slurm_cpus_per_task, stage_slurm_cfg.get("cpus_per_task")
        ),
        "slurm_output": pick(args.slurm_output, stage_slurm_cfg.get("output")),
        "slurm_error": pick(args.slurm_error, stage_slurm_cfg.get("error")),
        "submit": args.submit,
        "dry_run": args.dry_run,
        "_slurm_cfg_raw": slurm_cfg_raw,
    }

    settings["work_root"] = settings["work_root"] or "work"
    settings["fastp_threads"] = settings["fastp_threads"] or 8
    if settings["number_of_split_parts"] is not None:
        settings["number_of_split_parts"] = int(settings["number_of_split_parts"])
        if settings["number_of_split_parts"] <= 0:
            raise ValueError("number_of_split_parts must be > 0")
    settings["fastp_bin"] = normalize_executable_setting(
        settings["fastp_bin"], "fastp"
    )
    if stage == "demux_extract_bc":
        settings["barcode_whitelist"] = (
            settings["barcode_whitelist"] or DEFAULT_BARCODE_WHITELIST
        )
        settings["barcode_hamming_distance"] = int(
            settings["barcode_hamming_distance"] or 1
        )
        settings["gzip_level"] = int(settings["gzip_level"] or 6)
    settings["slurm_partition"] = settings["slurm_partition"] or "cpu"
    settings["slurm_mem"] = settings["slurm_mem"] or "16G"
    settings["slurm_cpus_per_task"] = settings["slurm_cpus_per_task"] or 8
    settings["slurm_output"] = settings["slurm_output"] or str(
        Path(settings["work_root"])
        / settings["sample_id"]
        / "logs"
        / f"{stage}_%x_%j.out"
    )
    settings["slurm_error"] = settings["slurm_error"] or str(
        Path(settings["work_root"])
        / settings["sample_id"]
        / "logs"
        / f"{stage}_%x_%j.err"
    )

    if stage == "all":
        for stage_name in STAGE_SEQUENCE:
            validate_required_for_stage(stage_name, settings)
    else:
        validate_required_for_stage(stage, settings)
    return settings


def generate_local_script(command: str, output_path: Path) -> None:
    content = (
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n\n"
        f"{command}\n"
    )
    write_text(output_path, content)
    output_path.chmod(0o755)


def generate_slurm_script(
    command: str, output_path: Path, log_dir: Path, args: argparse.Namespace
) -> None:
    log_dir.mkdir(parents=True, exist_ok=True)
    lines = [
        "#!/usr/bin/env bash",
        f"#SBATCH --job-name={args.job_name}",
        f"#SBATCH --partition={args.slurm_partition}",
        f"#SBATCH --cpus-per-task={args.slurm_cpus_per_task}",
        f"#SBATCH --mem={args.slurm_mem}",
        f"#SBATCH --output={args.slurm_output}",
        f"#SBATCH --error={args.slurm_error}",
        "",
        "set -euo pipefail",
        "",
        command,
        "",
    ]
    write_text(output_path, "\n".join(lines))
    output_path.chmod(0o755)


def submit_script(path: Path, runner: str) -> None:
    if runner == "local":
        subprocess.run(["bash", str(path)], check=True)
    else:
        subprocess.run(["sbatch", str(path)], check=True)


def parse_generated_paths(command_output: str) -> list[Path]:
    generated: list[Path] = []
    for line in command_output.splitlines():
        prefix = "[make_cmd] generated="
        if line.startswith(prefix):
            generated.append(Path(line[len(prefix) :].strip()))
    return generated


def build_stage_passthrough_args(argv: list[str]) -> list[str]:
    passthrough: list[str] = []
    flags_without_value = {
        "--submit",
        "--dry-run",
        "--skip-workdir-input-checks",
    }
    index = 0
    while index < len(argv):
        token = argv[index]
        if token in {"--stage", "--runner"}:
            index += 2
            continue
        if token.startswith("--stage=") or token.startswith("--runner="):
            index += 1
            continue
        if token in {"--submit", "--dry-run", "--skip-workdir-input-checks"}:
            index += 1
            continue
        if token in flags_without_value:
            passthrough.append(token)
            index += 1
            continue
        if token.startswith("--"):
            passthrough.append(token)
            if index + 1 < len(argv):
                passthrough.append(argv[index + 1])
            index += 2
            continue
        passthrough.append(token)
        index += 1
    return passthrough


def driver_scripts_for_stage(
    stage_name: str, scripts: list[Path], *, runner: str
) -> list[Path]:
    if stage_name != "demux_extract_bc":
        return scripts
    if runner == "local":
        return [script for script in scripts if script.name == "02_demux_extract_bc.sh"]
    return [script for script in scripts if script.suffix == ".sbatch"]


def generate_local_driver_script(
    stage_scripts: list[tuple[str, list[Path]]], output_path: Path
) -> None:
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        'SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
        "",
    ]
    for stage_name, scripts in stage_scripts:
        runnable = driver_scripts_for_stage(stage_name, scripts, runner="local")
        if not runnable:
            continue
        for script_path in runnable:
            lines.append(f'bash "$SCRIPT_DIR/{script_path.name}"')
    lines.append("")
    write_text(output_path, "\n".join(lines))
    output_path.chmod(0o755)


def generate_slurm_driver_script(
    stage_scripts: list[tuple[str, list[Path]]],
    output_path: Path,
    log_dir: Path,
    settings: dict,
) -> None:
    lines = [
        "submit_with_dep() {",
        '  local script_path="$1"',
        '  local dep_chain="$2"',
        "  local out",
        '  if [[ -n "$dep_chain" ]]; then',
        '    out="$(sbatch --dependency=afterok:${dep_chain} "$script_path")"',
        "  else",
        '    out="$(sbatch "$script_path")"',
        "  fi",
        '  echo "$out" >&2',
        '  echo "${out##* }"',
        "}",
        "",
        "join_deps() {",
        "  local joined=''",
        '  for item in "$@"; do',
        '    if [[ -z "$item" ]]; then',
        "      continue",
        "    fi",
        '    if [[ -z "$joined" ]]; then',
        '      joined="$item"',
        "    else",
        '      joined="${joined}:$item"',
        "    fi",
        "  done",
        '  echo "$joined"',
        "}",
        "",
        'SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
        'prev_stage_deps=""',
        "",
    ]
    for stage_name, scripts in stage_scripts:
        runnable = driver_scripts_for_stage(stage_name, scripts, runner="slurm")
        if not runnable:
            continue
        lines.append(f'echo "[run.sbatch] stage={stage_name}"')
        if stage_name == "demux_extract_bc":
            chunk_scripts = [
                script
                for script in runnable
                if script.name.startswith("02_demux_extract_bc_")
            ]
            aggregate_scripts = [
                script for script in runnable if script.name.startswith("02_aggregate_")
            ]
            chunk_job_vars: list[str] = []
            for index, script_path in enumerate(chunk_scripts):
                var_name = f"jid_demux_chunk_{index}"
                lines.append(
                    f'{var_name}="$(submit_with_dep "$SCRIPT_DIR/{script_path.name}" "$prev_stage_deps")"'
                )
                chunk_job_vars.append(var_name)
            if chunk_job_vars:
                deps_join = " ".join(f"${var_name}" for var_name in chunk_job_vars)
                lines.append(f'chunk_deps="$(join_deps {deps_join})"')
                for script_path in aggregate_scripts:
                    lines.append(
                        f'prev_stage_deps="$(submit_with_dep "$SCRIPT_DIR/{script_path.name}" "$chunk_deps")"'
                    )
            continue
        job_vars: list[str] = []
        for index, script_path in enumerate(runnable):
            var_name = f"jid_{stage_name}_{index}".replace("-", "_")
            lines.append(
                f'{var_name}="$(submit_with_dep "$SCRIPT_DIR/{script_path.name}" "$prev_stage_deps")"'
            )
            job_vars.append(var_name)
        deps_join = " ".join(f"${var_name}" for var_name in job_vars)
        lines.append(f'prev_stage_deps="$(join_deps {deps_join})"')
        lines.append("")
    lines.append('echo "[run.sbatch] done final_dep=${prev_stage_deps}"')
    driver_command = "\n".join(lines)
    slurm_args = argparse.Namespace(
        job_name=f"seeksoul_all_driver_{settings['sample_id']}",
        slurm_partition=settings["slurm_partition"],
        slurm_mem=settings["slurm_mem"],
        slurm_cpus_per_task=settings["slurm_cpus_per_task"],
        slurm_output=settings["slurm_output"].replace(
            "%x", f"seeksoul_all_driver_{settings['sample_id']}"
        ),
        slurm_error=settings["slurm_error"].replace(
            "%x", f"seeksoul_all_driver_{settings['sample_id']}"
        ),
    )
    generate_slurm_script(driver_command, output_path, log_dir, slurm_args)


def apply_stage_slurm_settings(settings: dict, stage: str) -> dict:
    slurm_cfg_raw = settings.get("_slurm_cfg_raw", {})
    stage_slurm_cfg = select_stage_slurm_cfg(slurm_cfg_raw, stage)
    updated = dict(settings)
    updated["slurm_partition"] = stage_slurm_cfg.get("partition") or settings["slurm_partition"]
    updated["slurm_mem"] = stage_slurm_cfg.get("mem") or settings["slurm_mem"]
    updated["slurm_cpus_per_task"] = (
        stage_slurm_cfg.get("cpus_per_task") or settings["slurm_cpus_per_task"]
    )
    updated["slurm_output"] = stage_slurm_cfg.get("output") or settings["slurm_output"]
    updated["slurm_error"] = stage_slurm_cfg.get("error") or settings["slurm_error"]
    return updated


def main() -> int:
    args = parse_args()
    settings = resolve_settings(args)
    sample_work = Path(settings["work_root"]) / settings["sample_id"]
    command_dir = sample_work / "commands"
    log_dir = sample_work / "logs"

    if settings["stage"] != "all":
        validate_inputs_for_stage(
            settings["stage"],
            settings,
            sample_work,
            skip_workdir_inputs=bool(args.skip_workdir_input_checks),
        )

    if settings["stage"] == "all":
        for stage_name in STAGE_SEQUENCE:
            validate_inputs_for_stage(
                stage_name,
                settings,
                sample_work,
                skip_workdir_inputs=True,
            )
        passthrough_args = build_stage_passthrough_args(sys.argv[1:])
        stage_scripts: list[tuple[str, list[Path]]] = []
        for stage_name in STAGE_SEQUENCE:
            stage_argv = [
                sys.executable,
                __file__,
                *passthrough_args,
                "--runner",
                settings["runner"],
                "--stage",
                stage_name,
            ]
            if settings["dry_run"]:
                stage_argv.append("--dry-run")
            stage_argv.append("--skip-workdir-input-checks")
            completed = subprocess.run(
                stage_argv, check=False, capture_output=True, text=True
            )
            if completed.stdout:
                print(completed.stdout, end="")
            if completed.stderr:
                print(completed.stderr, end="", file=sys.stderr)
            if completed.returncode != 0:
                return completed.returncode
            stage_scripts.append((stage_name, parse_generated_paths(completed.stdout)))

        driver_path: Path
        if settings["runner"] == "local":
            driver_path = command_dir / "run.sh"
            print(f"[make_cmd] script={driver_path}")
            if not settings["dry_run"]:
                generate_local_driver_script(stage_scripts, driver_path)
        else:
            driver_path = command_dir / "run.sbatch"
            print(f"[make_cmd] script={driver_path}")
            if not settings["dry_run"]:
                driver_settings = apply_stage_slurm_settings(
                    settings, STAGE_SEQUENCE[-1]
                )
                generate_slurm_driver_script(
                    stage_scripts, driver_path, log_dir, driver_settings
                )

        if not settings["dry_run"] and driver_path.exists():
            print(f"[make_cmd] generated={driver_path}")
        if settings["submit"] and not settings["dry_run"]:
            if settings["runner"] == "slurm":
                subprocess.run(["bash", str(driver_path)], check=True)
                print("[make_cmd] submitted_driver=1")
                print("[make_cmd] submit_mode=client_side_sbatch_dag")
            else:
                submit_script(driver_path, settings["runner"])
                print("[make_cmd] submitted_driver=1")

        print("[make_cmd] stage=all helper generation complete")
        return 0

    generated_scripts: list[Path] = []
    if settings["stage"] == "fastp_split":
        base_name = "01_fastp_split"
        command_args = argparse.Namespace(
            r1=settings["r1"],
            r2=settings["r2"],
            fastp_threads=settings["fastp_threads"],
            number_of_split_parts=settings["number_of_split_parts"],
            fastp_bin=settings["fastp_bin"],
        )
        command = build_fastp_split_command(command_args, sample_work)
        if settings["runner"] == "local":
            script_path = command_dir / f"{base_name}.sh"
        else:
            script_path = command_dir / f"{base_name}.sbatch"

        print(f"[make_cmd] runner={settings['runner']}")
        print(f"[make_cmd] stage={settings['stage']}")
        print(f"[make_cmd] sample_id={settings['sample_id']}")
        print(f"[make_cmd] script={script_path}")
        print(f"[make_cmd] command={command}")

        if settings["dry_run"]:
            return 0

        if settings["runner"] == "local":
            generate_local_script(command, script_path)
        else:
            slurm_args = argparse.Namespace(
                job_name=f"seeksoul_fastp_split_{settings['sample_id']}",
                slurm_partition=settings["slurm_partition"],
                slurm_mem=settings["slurm_mem"],
                slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                slurm_output=settings["slurm_output"],
                slurm_error=settings["slurm_error"],
            )
            generate_slurm_script(command, script_path, log_dir, slurm_args)
        generated_scripts.append(script_path)
    elif settings["stage"] == "demux_extract_bc":
        command_args = argparse.Namespace(
            barcode_whitelist=settings["barcode_whitelist"],
            barcode_hamming_distance=settings["barcode_hamming_distance"],
            gzip_level=settings["gzip_level"],
        )
        demux_dir = sample_work / "demux"
        if settings["runner"] == "local":
            script_path = command_dir / "02_demux_extract_bc.sh"
            command = build_demux_local_batch_command(
                command_args, sample_work, settings.get("number_of_split_parts")
            )
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] script={script_path}")
            print(f"[make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            generated_scripts.append(script_path)
        else:
            chunks = discover_demux_chunks(
                sample_work, settings.get("number_of_split_parts")
            )
            if not chunks:
                raise ValueError("no fastp shards found for demux script generation")
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] chunk_count={len(chunks)}")
            for chunk_id, r1_path, r2_path, out_prefix in chunks:
                base_name = f"02_demux_extract_bc_{chunk_id}"
                script_path = command_dir / f"{base_name}.sbatch"
                command = build_demux_chunk_command(
                    command_args, r1_path, r2_path, out_prefix
                )
                chunk_output = settings["slurm_output"].replace(
                    "%x", f"seeksoul_demux_{settings['sample_id']}_{chunk_id}"
                )
                chunk_error = settings["slurm_error"].replace(
                    "%x", f"seeksoul_demux_{settings['sample_id']}_{chunk_id}"
                )
                print(f"[make_cmd] script={script_path}")
                print(f"[make_cmd] command={command}")
                if not settings["dry_run"]:
                    slurm_args = argparse.Namespace(
                        job_name=f"seeksoul_demux_{settings['sample_id']}_{chunk_id}",
                        slurm_partition=settings["slurm_partition"],
                        slurm_mem=settings["slurm_mem"],
                        slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                        slurm_output=chunk_output,
                        slurm_error=chunk_error,
                    )
                    generate_slurm_script(command, script_path, log_dir, slurm_args)
                generated_scripts.append(script_path)

            aggregate_script = command_dir / "02_aggregate_ct_qc.sbatch"
            aggregate_command = build_aggregate_ct_command(demux_dir)
            aggregate_output = settings["slurm_output"].replace(
                "%x", f"seeksoul_aggregate_ct_{settings['sample_id']}"
            )
            aggregate_error = settings["slurm_error"].replace(
                "%x", f"seeksoul_aggregate_ct_{settings['sample_id']}"
            )
            print(f"[make_cmd] script={aggregate_script}")
            print(f"[make_cmd] command={aggregate_command}")
            if not settings["dry_run"]:
                slurm_args = argparse.Namespace(
                    job_name=f"seeksoul_aggregate_ct_{settings['sample_id']}",
                    slurm_partition=settings["slurm_partition"],
                    slurm_mem=settings["slurm_mem"],
                    slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                    slurm_output=aggregate_output,
                    slurm_error=aggregate_error,
                )
                generate_slurm_script(
                    aggregate_command, aggregate_script, log_dir, slurm_args
                )
                submit_script_path = command_dir / "02_demux_extract_bc_submit.sh"
                generate_local_script(
                    build_demux_slurm_submit_command(sample_work),
                    submit_script_path,
                )
                generated_scripts.append(aggregate_script)
                generated_scripts.append(submit_script_path)
            else:
                generated_scripts.append(aggregate_script)
    else:
        raise ValueError(f"unsupported stage: {settings['stage']}")

    for script_path in generated_scripts:
        print(f"[make_cmd] generated={script_path}")

    if settings["submit"]:
        for script_path in generated_scripts:
            submit_script(script_path, settings["runner"])
        print(f"[make_cmd] submitted_count={len(generated_scripts)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
