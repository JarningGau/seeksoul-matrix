## HPC test
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dd_met5_test.json \
  --stage all \
  --runner slurm \
  --r1 /storage2/liliLab/gaojianing/SeekSpace-DNAme-data/test_R1.fastq.gz \
  --r2 /storage2/liliLab/gaojianing/SeekSpace-DNAme-data/test_R2.fastq.gz \
  --sample-id test \
  --bismark-ref /storage2/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/fasta/ \
  --genome-fa /storage2/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/fasta/genome.fa \
  --chrom-size-path /storage2/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/bed/chr_nochrM.bed \
  --submit

## HPC C283
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dd_met5_slurm.json \
  --stage all \
  --split-fastq-prefix-bases 2 \
  --force-cell-num 10000 \
  --r1 /storage2/liliLab/gaojianing/SeekSpace-DNAme-data/C283_Brain_DNAme_R1.fastq.gz \
  --r2 /storage2/liliLab/gaojianing/SeekSpace-DNAme-data/C283_Brain_DNAme_R2.fastq.gz \
  --sample-id C283_Brain_DNAme \
  --bismark-ref /storage2/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/fasta/ \
  --genome-fa /storage2/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/fasta/genome.fa \
  --chrom-size-path /storage2/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/bed/chr_nochrM.bed


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
--stage meth_matrix \
--sample-id C283_Brain_DNAme \
--meth-regions-bed /storage2/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/regions/bin500.bed \
--meth-regions-label bin500

## HPC C283 multi-lane (copy examples/lanes.tsv.example -> examples/lanes.tsv)
## Login node: submit one uploaded lane through bismark, wait for jobs, harvest, repeat.
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dd_met5_slurm.json \
  --stage fastp_split demux_extract_bc regroup_shards bismark_align \
  --runner slurm \
  --sample-id L01 \
  --work-root work/C283_Brain_DNAme_lanes \
  --r1 /storage2/liliLab/gaojianing/SeekSpace-DNAme-data/C283_Brain/C283_Brain_DNAme_pilot_R1.fastq.gz \
  --r2 /storage2/liliLab/gaojianing/SeekSpace-DNAme-data/C283_Brain/C283_Brain_DNAme_pilot_R2.fastq.gz \
  --bismark-ref /storage2/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/fasta/ \
  --genome-fa /storage2/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/fasta/genome.fa \
  --chrom-size-path /storage2/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/bed/chr_nochrM.bed \
  --submit

bash examples/run_multi_lane.sh \
  --manifest examples/lanes.tsv \
  --sample-id C283_Brain_DNAme \
  --phase harvest

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dd_met5_slurm.json \
  --stage fastp_split demux_extract_bc regroup_shards bismark_align \
  --runner slurm \
  --sample-id L02 \
  --work-root work/C283_Brain_DNAme_lanes \
  --r1 /storage2/liliLab/gaojianing/SeekSpace-DNAme-data/C283_Brain/202609121633_B_SG260912-6_XYRD_WTJW1090_Met_L00_R1.fq.gz \
  --r2 /storage2/liliLab/gaojianing/SeekSpace-DNAme-data/C283_Brain/202609121633_B_SG260912-6_XYRD_WTJW1090_Met_L00_R2.fq.gz \
  --bismark-ref /storage2/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/fasta/ \
  --genome-fa /storage2/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/fasta/genome.fa \
  --chrom-size-path /storage2/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/bed/chr_nochrM.bed \
  --submit

bash examples/run_multi_lane.sh \
  --manifest examples/lanes.tsv \
  --sample-id C283_Brain_DNAme \
  --phase harvest

## After every manifest lane has work/C283_Brain_DNAme/lanes/<id>.done
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dd_met5_slurm.json \
  --stage bam_sort count_mapped_reads estimated_cells split_bams merge_fr_bams bam_to_allc saturation qc_summary \
  --runner slurm \
  --sample-id C283_Brain_DNAme \
  --work-root work \
  --force-cell-num 10000 \
  --bismark-ref /storage2/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/fasta/ \
  --genome-fa /storage2/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/fasta/genome.fa \
  --chrom-size-path /storage2/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/bed/chr_nochrM.bed \
  --submit
