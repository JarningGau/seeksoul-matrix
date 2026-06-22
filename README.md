# SeepSpace Spatial Methylation (Slide-Tag like)

## DD-MET5 Chemistry

### Library structure
`CB`-`UB`-`TSO`-`17L`-`ME`-`9bp`-`Insert`-`9bp`

- `CB`: 17 bp cell barcode (do not contain C)
- `UB`: 12 bp UMI sequence (do not contain C)
- `TSO`: 13 bp TSO sequence TTTCTTATATGGG
- `17L`: 17 bp fixed sequence CgtCCgtCgttgCtCgt
- `ME`: 19 bp fixed sequence AGATGTGTATAAGAGACAG
- `9bp`: extension sequence from the Tn5 insertion fragment

Since the enzymatic treatment converts unmethylated cytosines (C) to thymines (T), the barcodes used for methylation data do not contain any C bases. In contrast, the C bases in TSO, 17L, and ME are not methylated and will be converted to T during the enzymatic process; we use these fixed sequences to calculate the C-to-T conversion rate.

### Cell Barcode Whitelist Validation

**Whitelist**: `whitelist/DD-MET5/U3CB_methylation.txt.gz` (829,440 × 17 bp cell barcodes)

### Installation

Requires [pixi](https://pixi.sh/) on **linux-64**.

```bash
git clone https://github.com/JarningGau/seeksoul-matrix.git
cd seeksoul-matrix
pixi install
```

The pixi environment provides core tools (`fastp`, `cutadapt`, `bowtie2`, `samtools`, Python 3.11). Two seekgene forks are installed separately into the active environment (not conda packages):

```bash
# bismark_align — seekgene/Bismark with --add_barcode / --add_umi
pixi run setup-bismark
pixi run check-bismark-env

# bam_to_allc — seekgene/ALLCools with UR-tag UMI dedup
pixi run setup-allcools
pixi run check-allcools-env
```

Re-run `setup-bismark` and `setup-allcools` after `pixi install` recreates `.pixi/envs/default`.
