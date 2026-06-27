## HPC test
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dd_met5_test.json \
  --stage all \
  --runner slurm \
  --r1 /storage/liliLab/gaojianing/SeekSpace-DNAme-data/test_R1.fastq.gz \
  --r2 /storage/liliLab/gaojianing/SeekSpace-DNAme-data/test_R2.fastq.gz \
  --sample-id test \
  --bismark-ref /storage/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/fasta/ \
  --genome-fa /storage/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/fasta/genome.fa \
  --chrom-size-path /storage/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/bed/chr_nochrM.bed \
  --submit

## HPC C283
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dd_met5_slurm.json \
  --stage all \
  --split-fastq-prefix-bases 2 \
  --force-cell-num 10000 \
  --r1 /storage/liliLab/gaojianing/SeekSpace-DNAme-data/C283_Brain_DNAme_R1.fastq.gz \
  --r2 /storage/liliLab/gaojianing/SeekSpace-DNAme-data/C283_Brain_DNAme_R2.fastq.gz \
  --sample-id C283_Brain_DNAme \
  --bismark-ref /storage/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/fasta/ \
  --genome-fa /storage/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/fasta/genome.fa \
  --chrom-size-path /storage/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/bed/chr_nochrM.bed


## HPC C283 qc_summary only
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dd_met5_slurm.json \
  --stage qc_summary \
  --sample-id C283_Brain_DNAme

## HPC C283 meth analysis
for stage in allc_to_matrix meth_smooth meth_scan meth_matrix; do
  pixi run python scripts/make_cmd.py \
    --workflow-config workflow/dd_met5_slurm.json \
    --stage "$stage" \
    --sample-id C283_Brain_DNAme
done

pixi run python scripts/make_cmd.py \
--workflow-config workflow/dd_met5_slurm.json \
--stage meth_smooth \
--sample-id C283_Brain_DNAme
