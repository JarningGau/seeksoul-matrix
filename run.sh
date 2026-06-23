pixi run python scripts/make_cmd.py --workflow-config workflow/dd_met5_test.json --stage all \
--r1 /mnt/e/Spatio_DARLIN_data/SeekSpace/DNAme/fastq/C283_Brain_DNAme_S1_R1.fastq.gz \
--r2 /mnt/e/Spatio_DARLIN_data/SeekSpace/DNAme/fastq/C283_Brain_DNAme_S1_R2.fastq.gz \
--sample-id C283_Brain_DNAme_S1 --submit


pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dd_met5_test.json \
  --stage all \
  --runner slurm \
  --number-of-split-parts 2 \
  --r1 /storage/liliLab/gaojianing/SeekSpace-DNAme-data/test_R1.fastq.gz \
  --r2 /storage/liliLab/gaojianing/SeekSpace-DNAme-data/test_R2.fastq.gz \
  --sample-id test \
  --bismark-ref /storage/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/fasta/ \
  --genome-fa /storage/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/fasta/genome.fa \
  --chrom-size-path /storage/liliLab/gaojianing/resource/seeksoul/mouse-reference-GRCm39/bed/chr_nochrM.bed \
  --submit