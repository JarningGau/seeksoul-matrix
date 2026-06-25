## Local test
pixi run python scripts/make_cmd.py --workflow-config workflow/dd_met5_test.json --stage all \
--r1 /mnt/e/Spatio_DARLIN_data/SeekSpace/DNAme/fastq/C283_Brain_DNAme_S1_R1.fastq.gz \
--r2 /mnt/e/Spatio_DARLIN_data/SeekSpace/DNAme/fastq/C283_Brain_DNAme_S1_R2.fastq.gz \
--sample-id C283_Brain_DNAme_S1 --submit
