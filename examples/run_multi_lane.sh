#!/usr/bin/env bash
# Disk-aware multi-lane queue around the existing seeksoul-matrix stages.
#
# Does not change stage scripts or make_cmd.py contracts. Each lane is a
# temporary sample through bismark_align; unsorted BAMs are concatenated into
# one final sample; later stages run there with explicit rm between copies.
#
# HPC alignment: on the login node, submit each uploaded lane with
# make_cmd --runner slurm --submit and a contiguous --stage list through
# bismark_align (see examples/run_slurm_example.sh). After that lane's jobs
# finish, run this wrapper with --phase harvest.
#
# Local / non-Slurm: --phase lanes|all still runs make_cmd --runner local.
#
# Usage (from repository root):
#   bash examples/run_multi_lane.sh \
#     --manifest examples/lanes.tsv \
#     --sample-id C283_Brain_DNAme \
#     --workflow-config workflow/dd_met5_slurm.json
#
# Rolling upload (one lane on disk at a time):
#   bash examples/run_multi_lane.sh \
#     --manifest examples/lanes.tsv \
#     --sample-id C283_Brain_DNAme \
#     --skip-missing \
#     --delete-raw-fastq
# Re-run the same command after each upload. Tail stages (cell calling / ALLC)
# start only when every manifest lane has work/<final>/lanes/<id>.done.
#
# Resume: skip lanes with work/<final>/lanes/<lane_id>.done
# If align/.merging exists, the previous BAM cat crashed — see the file
# and restore or delete the destination BAM before retrying.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

usage() {
  cat <<'EOF'
Usage: bash examples/run_multi_lane.sh --manifest PATH --sample-id ID [options]

Required:
  --manifest PATH           TSV of lane_id, r1, r2 (see examples/lanes.tsv.example)
  --sample-id ID            Final sample directory under --work-root

Options:
  --workflow-config PATH    Default: workflow/dd_met5_slurm.json
  --work-root PATH          Default: work
  --phase all|lanes|harvest|qc|tail Default: all
  --dry-run                 Print make_cmd / rm / samtools plans only
  --skip-missing            Skip lanes whose FASTQs are not on disk yet
  --delete-raw-fastq        Delete source R1/R2 after that lane is merged
  --keep-lane-work          Keep per-lane work dirs after merge
  --keep-unsorted-bam       Keep unsorted Bismark BAMs after bam_sort
  --keep-sortbyname-bam     Keep name-sorted BAMs after split_bams
  --keep-split-strand-bams  Keep per-strand split BAM dirs after merge_fr
  --keep-merged-cell-bams   Keep merged per-cell BAMs after saturation
  --samtools-bin PATH       samtools executable (default: samtools or pixi run)
  -h, --help                Show this help

Submit --phase harvest after a lane's Slurm DAG finishes. Local --phase lanes|all still use make_cmd --runner local.
EOF
}

need_arg() {
  if [[ $# -lt 2 || -z "${2}" || "$2" == --* ]]; then
    echo "[run_multi_lane] missing value for $1" >&2
    exit 2
  fi
}

LANE_MANIFEST=""
FINAL_SAMPLE=""
WORKFLOW_CONFIG="workflow/dd_met5_slurm.json"
WORK_ROOT="work"
PHASE="all"
DRY_RUN=0
SKIP_MISSING=0
DELETE_LANE_WORK=1
DELETE_RAW_FASTQ=0
DELETE_UNSORTED_BAM=1
DELETE_SORTBYNAME_BAM=1
DELETE_SPLIT_STRAND_BAMS=1
DELETE_MERGED_CELL_BAMS=1
SAMTOOLS_BIN=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest)
      need_arg "$@"
      LANE_MANIFEST="$2"
      shift 2
      ;;
    --sample-id)
      need_arg "$@"
      FINAL_SAMPLE="$2"
      shift 2
      ;;
    --workflow-config)
      need_arg "$@"
      WORKFLOW_CONFIG="$2"
      shift 2
      ;;
    --work-root)
      need_arg "$@"
      WORK_ROOT="$2"
      shift 2
      ;;
    --phase)
      need_arg "$@"
      PHASE="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    --skip-missing)
      SKIP_MISSING=1
      shift
      ;;
    --delete-raw-fastq)
      DELETE_RAW_FASTQ=1
      shift
      ;;
    --keep-lane-work)
      DELETE_LANE_WORK=0
      shift
      ;;
    --keep-unsorted-bam)
      DELETE_UNSORTED_BAM=0
      shift
      ;;
    --keep-sortbyname-bam)
      DELETE_SORTBYNAME_BAM=0
      shift
      ;;
    --keep-split-strand-bams)
      DELETE_SPLIT_STRAND_BAMS=0
      shift
      ;;
    --keep-merged-cell-bams)
      DELETE_MERGED_CELL_BAMS=0
      shift
      ;;
    --samtools-bin)
      need_arg "$@"
      SAMTOOLS_BIN="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      break
      ;;
    -*)
      echo "[run_multi_lane] unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
    *)
      echo "[run_multi_lane] unexpected argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -z "$LANE_MANIFEST" || -z "$FINAL_SAMPLE" ]]; then
  echo "[run_multi_lane] --manifest and --sample-id are required" >&2
  usage >&2
  exit 2
fi

MAKE_CMD=(pixi run python scripts/make_cmd.py)
SAMTOOLS=(samtools)
if [[ -n "$SAMTOOLS_BIN" ]]; then
  SAMTOOLS=("$SAMTOOLS_BIN")
elif ! command -v samtools >/dev/null 2>&1; then
  SAMTOOLS=(pixi run samtools)
fi

FINAL_WORK="$WORK_ROOT/$FINAL_SAMPLE"
LANE_WORK_ROOT="$WORK_ROOT/${FINAL_SAMPLE}_lanes"
FINAL_ALIGN="$FINAL_WORK/align"
FINAL_DEMUX="$FINAL_WORK/demux"
LANE_QC_ROOT="$FINAL_WORK/lane_qc"
LANE_DONE_DIR="$FINAL_WORK/lanes"
MERGE_LOCK="$FINAL_ALIGN/.merging"

if [[ "$DRY_RUN" != "1" ]]; then
  mkdir -p "$FINAL_WORK" "$LANE_WORK_ROOT" "$FINAL_ALIGN" "$FINAL_DEMUX" \
    "$LANE_QC_ROOT" "$LANE_DONE_DIR" "$FINAL_WORK/shard_fastq"
fi

log() { echo "[run_multi_lane] $*"; }

run() {
  if [[ "$DRY_RUN" == "1" ]]; then
    log "dry-run: $*"
    return 0
  fi
  "$@"
}

make_cmd_stage() {
  local sample_id="$1"
  local stage="$2"
  shift 2
  local extra=("$@")
  local cmd=(
    "${MAKE_CMD[@]}"
    --workflow-config "$WORKFLOW_CONFIG"
    --runner local
    --sample-id "$sample_id"
    --stage "$stage"
  )
  cmd+=("${extra[@]}")
  if [[ "$DRY_RUN" == "1" ]]; then
    cmd+=(--dry-run --skip-workdir-input-checks)
    log "dry-run make_cmd sample=$sample_id stage=$stage"
    "${cmd[@]}"
    return 0
  fi
  cmd+=(--submit)
  log "make_cmd sample=$sample_id stage=$stage"
  "${cmd[@]}"
}

fastq_ready() {
  local r1="$1" r2="$2"
  [[ -f "$r1" && -s "$r1" && -f "$r2" && -s "$r2" ]]
}

lane_is_done() {
  [[ -f "$LANE_DONE_DIR/${1}.done" ]]
}

count_done_lanes() {
  local lane_id n=0
  for lane_id in "${LANES[@]}"; do
    if lane_is_done "$lane_id"; then
      n=$((n + 1))
    fi
  done
  echo "$n"
}

valid_lane_id() {
  [[ "$1" =~ ^[A-Za-z0-9-]+$ ]]
}

read_manifest() {
  local line lane_id r1 r2
  LANES=()
  R1_PATHS=()
  R2_PATHS=()
  while IFS= read -r line || [[ -n "$line" ]]; do
    line="${line%%$'\r'}"
    [[ -z "$line" || "$line" == \#* ]] && continue
    IFS=$'\t' read -r lane_id r1 r2 _ <<<"$line"
    if [[ -z "$lane_id" || -z "$r1" || -z "$r2" ]]; then
      echo "[run_multi_lane] bad manifest row: $line" >&2
      exit 1
    fi
    if ! valid_lane_id "$lane_id"; then
      echo "[run_multi_lane] lane_id must match [A-Za-z0-9-]+: $lane_id" >&2
      exit 1
    fi
    LANES+=("$lane_id")
    R1_PATHS+=("$r1")
    R2_PATHS+=("$r2")
  done < "$LANE_MANIFEST"
  if [[ ${#LANES[@]} -eq 0 ]]; then
    echo "[run_multi_lane] no lanes in $LANE_MANIFEST" >&2
    exit 1
  fi
}

lane_work() { echo "$LANE_WORK_ROOT/$1"; }

has_unsorted_bams() {
  local dir="$1"
  local match
  shopt -s nullglob
  match=("$dir"/*.forward_1_bismark_bt2_pe.bam)
  shopt -u nullglob
  [[ ${#match[@]} -gt 0 ]]
}

run_lane_stages() {
  local lane_id="$1" r1="$2" r2="$3"
  local extra=(--work-root "$LANE_WORK_ROOT" --r1 "$r1" --r2 "$r2")
  local work
  work="$(lane_work "$lane_id")"
  if has_unsorted_bams "$work/align"; then
    log "lane=$lane_id skip stages (align BAMs present)"
    return 0
  fi
  make_cmd_stage "$lane_id" fastp_split "${extra[@]}"
  make_cmd_stage "$lane_id" demux_extract_bc "${extra[@]}"
  make_cmd_stage "$lane_id" regroup_shards "${extra[@]}"
  make_cmd_stage "$lane_id" bismark_align "${extra[@]}"
}

harvest_lane_qc() {
  local lane_id="$1"
  local work dest
  work="$(lane_work "$lane_id")"
  dest="$LANE_QC_ROOT/$lane_id"
  run mkdir -p "$dest/demux" "$dest/align" "$dest/shard_fastq"
  if [[ "$DRY_RUN" == "1" ]]; then
    log "dry-run harvest qc lane=$lane_id"
    return 0
  fi
  if [[ -f "$work/shard_fastq/fastp.json" ]]; then
    cp -a "$work/shard_fastq/fastp.json" "$dest/shard_fastq/fastp.json"
  fi
  local src
  shopt -s nullglob
  for src in "$work/demux"/*.linker.tsv "$work/demux"/*.stats.json; do
    cp -a "$src" "$dest/demux/${lane_id}_$(basename "$src")"
  done
  for src in "$work/align"/*_bismark_bt2_PE_report.txt; do
    cp -a "$src" "$dest/align/${lane_id}_$(basename "$src")"
  done
  shopt -u nullglob
}

merge_lane_bams() {
  local lane_id="$1"
  local work bam dest tmp name
  work="$(lane_work "$lane_id")"
  if [[ -f "$MERGE_LOCK" ]]; then
    echo "[run_multi_lane] incomplete BAM cat: $MERGE_LOCK" >&2
    echo "[run_multi_lane] inspect the lock, fix $FINAL_ALIGN, then remove the lock." >&2
    exit 1
  fi
  if [[ "$DRY_RUN" == "1" ]]; then
    log "dry-run samtools cat lane=$lane_id -> $FINAL_ALIGN"
    return 0
  fi
  if ! has_unsorted_bams "$work/align"; then
    echo "[run_multi_lane] no Bismark PE BAMs in $work/align" >&2
    exit 1
  fi
  mkdir -p "$FINAL_ALIGN"
  printf '%s\n' "$lane_id" > "$MERGE_LOCK"
  shopt -s nullglob
  for bam in "$work/align"/*.forward_1_bismark_bt2_pe.bam \
             "$work/align"/*.reverse_1_bismark_bt2_pe.bam; do
    name="$(basename "$bam")"
    dest="$FINAL_ALIGN/$name"
    tmp="$FINAL_ALIGN/${name}.partial"
    rm -f "$tmp"
    if [[ -f "$dest" ]]; then
      log "samtools cat $dest + $name"
      "${SAMTOOLS[@]}" cat -o "$tmp" "$dest" "$bam"
      mv -f "$tmp" "$dest"
    else
      log "mv $name -> $FINAL_ALIGN"
      mv -f "$bam" "$dest"
      continue
    fi
    rm -f "$bam"
  done
  shopt -u nullglob
  rm -f "$MERGE_LOCK"
}

cleanup_lane() {
  local lane_id="$1" r1="$2" r2="$3"
  local work
  work="$(lane_work "$lane_id")"
  if [[ "$DELETE_LANE_WORK" == "1" ]]; then
    log "rm lane work $work"
    run rm -rf "$work"
  else
    log "rm FASTQ intermediates in $work"
    if [[ "$DRY_RUN" == "1" ]]; then
      log "dry-run rm FASTQ intermediates in $work"
    else
      find "$work/shard_fastq" "$work/demux" -type f \( -name '*.fq.gz' -o -name '*.fastq.gz' \) -delete 2>/dev/null || true
      rm -rf "$work/demux/shards"
    fi
  fi
  if [[ "$DELETE_RAW_FASTQ" == "1" ]]; then
    log "rm raw FASTQ $r1 $r2"
    run rm -f "$r1" "$r2"
  fi
}

mark_lane_done() {
  local lane_id="$1"
  if [[ "$DRY_RUN" == "1" ]]; then
    log "dry-run mark done $lane_id"
    return 0
  fi
  date -Iseconds > "$LANE_DONE_DIR/${lane_id}.done"
}

mark_lane_aligned() {
  local lane_id="$1"
  if [[ "$DRY_RUN" == "1" ]]; then
    log "dry-run mark aligned $lane_id"
    return 0
  fi
  date -Iseconds > "$LANE_DONE_DIR/${lane_id}.aligned"
}

process_lanes() {
  local i lane_id r1 r2 work
  for i in "${!LANES[@]}"; do
    lane_id="${LANES[$i]}"
    r1="${R1_PATHS[$i]}"
    r2="${R2_PATHS[$i]}"
    work="$(lane_work "$lane_id")"
    if lane_is_done "$lane_id"; then
      log "lane=$lane_id already merged; skip"
      if [[ "$DELETE_LANE_WORK" == "1" && -d "$work" ]]; then
        cleanup_lane "$lane_id" "$r1" "$r2"
      fi
      continue
    fi
    if ! fastq_ready "$r1" "$r2"; then
      if [[ "$DRY_RUN" == "1" ]]; then
        log "lane=$lane_id FASTQ not on disk (dry-run continues)"
      elif [[ "$SKIP_MISSING" == "1" ]]; then
        log "lane=$lane_id FASTQ not on disk yet; skip"
        continue
      else
        echo "[run_multi_lane] missing FASTQ for lane=$lane_id" >&2
        echo "[run_multi_lane] r1=$r1 r2=$r2" >&2
        echo "[run_multi_lane] upload this pair first, or pass --skip-missing" >&2
        exit 1
      fi
    fi
    if [[ -f "$LANE_DONE_DIR/${lane_id}.aligned" ]] \
      && ! has_unsorted_bams "$work/align" \
      && [[ "$DRY_RUN" != "1" ]]; then
      echo "[run_multi_lane] lane=$lane_id aligned but BAMs are gone and .done is missing." >&2
      echo "[run_multi_lane] BAM cat likely crashed after updating $FINAL_ALIGN; do not re-run this lane." >&2
      exit 1
    fi
    log "lane=$lane_id start r1=$r1"
    run_lane_stages "$lane_id" "$r1" "$r2"
    mark_lane_aligned "$lane_id"
    harvest_lane_qc "$lane_id"
    merge_lane_bams "$lane_id"
    mark_lane_done "$lane_id"
    cleanup_lane "$lane_id" "$r1" "$r2"
    log "lane=$lane_id done"
  done
}

harvest_finished_lanes() {
  local i lane_id r1 r2 work
  for i in "${!LANES[@]}"; do
    lane_id="${LANES[$i]}"
    r1="${R1_PATHS[$i]}"
    r2="${R2_PATHS[$i]}"
    work="$(lane_work "$lane_id")"
    if lane_is_done "$lane_id"; then
      log "lane=$lane_id already merged; skip"
      if [[ "$DELETE_LANE_WORK" == "1" && -d "$work" ]]; then
        cleanup_lane "$lane_id" "$r1" "$r2"
      fi
      continue
    fi
    if ! has_unsorted_bams "$work/align"; then
      log "lane=$lane_id no unsorted BAMs yet; skip"
      continue
    fi
    log "lane=$lane_id harvest r1=$r1"
    mark_lane_aligned "$lane_id"
    harvest_lane_qc "$lane_id"
    merge_lane_bams "$lane_id"
    mark_lane_done "$lane_id"
    cleanup_lane "$lane_id" "$r1" "$r2"
    log "lane=$lane_id done"
  done
}

finalize_qc() {
  local lane_id src dest
  log "assemble sample-level QC glue"
  if [[ "$DRY_RUN" == "1" ]]; then
    log "dry-run finalize qc"
    return 0
  fi
  mkdir -p "$FINAL_DEMUX" "$FINAL_ALIGN" "$FINAL_WORK/shard_fastq"
  shopt -s nullglob
  for src in "$LANE_QC_ROOT"/*/demux/*.linker.tsv \
             "$LANE_QC_ROOT"/*/demux/*.stats.json; do
    dest="$FINAL_DEMUX/$(basename "$src")"
    cp -a "$src" "$dest"
  done
  for src in "$LANE_QC_ROOT"/*/align/*_bismark_bt2_PE_report.txt; do
    cp -a "$src" "$FINAL_ALIGN/$(basename "$src")"
  done
  local fastp_inputs=("$LANE_QC_ROOT"/*/shard_fastq/fastp.json)
  shopt -u nullglob
  if [[ ${#fastp_inputs[@]} -eq 0 ]]; then
    echo "[run_multi_lane] no lane fastp.json under $LANE_QC_ROOT" >&2
    exit 1
  fi
  pixi run python examples/merge_fastp_json.py \
    "${fastp_inputs[@]}" \
    -o "$FINAL_WORK/shard_fastq/fastp.json"
  pixi run python scripts/aggregate_ct_qc.py --demux-dir "$FINAL_DEMUX"
}

run_tail() {
  local extra=(--work-root "$WORK_ROOT")
  make_cmd_stage "$FINAL_SAMPLE" count_mapped_reads "${extra[@]}"
  make_cmd_stage "$FINAL_SAMPLE" estimated_cells "${extra[@]}"
  make_cmd_stage "$FINAL_SAMPLE" bam_sort "${extra[@]}"
  if [[ "$DELETE_UNSORTED_BAM" == "1" ]]; then
    log "rm unsorted Bismark BAMs"
    if [[ "$DRY_RUN" != "1" ]]; then
      find "$FINAL_ALIGN" -maxdepth 1 -type f -name '*_bismark_bt2_pe.bam' \
        ! -name '*sortbyname*' -delete
    fi
  fi
  make_cmd_stage "$FINAL_SAMPLE" split_bams "${extra[@]}"
  if [[ "$DELETE_SORTBYNAME_BAM" == "1" ]]; then
    log "rm sortbyname BAMs"
    if [[ "$DRY_RUN" == "1" ]]; then
      log "dry-run rm $FINAL_ALIGN/*_sortbyname.bam"
    else
      find "$FINAL_ALIGN" -maxdepth 1 -type f -name '*_sortbyname.bam' -delete
    fi
  fi
  make_cmd_stage "$FINAL_SAMPLE" merge_fr_bams "${extra[@]}"
  if [[ "$DELETE_SPLIT_STRAND_BAMS" == "1" ]]; then
    log "rm split strand BAM dirs"
    if [[ "$DRY_RUN" == "1" ]]; then
      log "dry-run rm $FINAL_WORK/split_bams/*.{forward,reverse}_1"
    else
      shopt -s nullglob
      local strand_dir
      for strand_dir in "$FINAL_WORK/split_bams"/*.forward_1 \
                        "$FINAL_WORK/split_bams"/*.reverse_1; do
        rm -rf "$strand_dir"
      done
      shopt -u nullglob
    fi
  fi
  make_cmd_stage "$FINAL_SAMPLE" bam_to_allc "${extra[@]}"
  make_cmd_stage "$FINAL_SAMPLE" saturation "${extra[@]}"
  if [[ "$DELETE_MERGED_CELL_BAMS" == "1" ]]; then
    log "rm merged per-cell BAMs"
    run rm -rf "$FINAL_WORK/split_bams/merged"
  fi
  make_cmd_stage "$FINAL_SAMPLE" qc_summary "${extra[@]}"
}

read_manifest
log "manifest=$LANE_MANIFEST lanes=${#LANES[@]} final=$FINAL_SAMPLE phase=$PHASE dry_run=$DRY_RUN skip_missing=$SKIP_MISSING"

case "$PHASE" in
  lanes)
    process_lanes
    if [[ "$(count_done_lanes)" -gt 0 ]]; then
      finalize_qc
    else
      log "no merged lanes yet; skip QC glue"
    fi
    ;;
  harvest)
    harvest_finished_lanes
    if [[ "$(count_done_lanes)" -gt 0 ]]; then
      finalize_qc
    else
      log "no merged lanes yet; skip QC glue"
    fi
    ;;
  qc)
    finalize_qc
    ;;
  tail)
    if [[ "$(count_done_lanes)" -ne "${#LANES[@]}" ]]; then
      echo "[run_multi_lane] tail needs every manifest lane merged ($(count_done_lanes)/${#LANES[@]} done)" >&2
      exit 1
    fi
    run_tail
    ;;
  all)
    process_lanes
    if [[ "$(count_done_lanes)" -gt 0 ]]; then
      finalize_qc
    else
      log "no merged lanes yet; skip QC glue"
    fi
    if [[ "$(count_done_lanes)" -eq "${#LANES[@]}" ]]; then
      run_tail
    else
      log "merged $(count_done_lanes)/${#LANES[@]} lanes; skip tail until the rest are uploaded"
    fi
    ;;
  *)
    echo "[run_multi_lane] --phase must be all|lanes|harvest|qc|tail (got $PHASE)" >&2
    exit 1
    ;;
esac

log "phase=$PHASE complete"
