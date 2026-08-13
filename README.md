# xpclr-cpp

<!-- README-I18N:START -->

**English** | [汉语](./README.zh.md)

<!-- README-I18N:END -->

C++ XP-CLR (Chen, Patterson & Reich 2010) with htslib VCF/BCF I/O.

Compared with [hardingnj/xpclr](https://github.com/hardingnj/xpclr):

- region queries via TBI/CSI (`-r`)
- one `SAMPLE GROUP` pop file (`-a` / `-b`)
- full selection-coefficient grid by default (`--unimodal-s` for Python early-stop)
- `--seed` for `maxsnps` subsample; `--omega-trim` for background ω; optional `--gmap`
- `--phased INT` picks how LD weights handle missing genotypes and phasing

## Install

```bash
git clone --recurse-submodules https://github.com/WWz33/xpclr-cpp.git
cd xpclr-cpp
make -j
```

Needs a C++17 compiler, zlib, bzip2, lzma, libcurl, openssl, libdeflate.
Vendored htslib and GSL build with `make`. System libs: `make USE_SYSTEM_HTS=1 USE_SYSTEM_GSL=1`.

## Example

```bash
./xpclr -i data/smoke.vcf.gz -p data/pop_smoke.txt \
  -a popA -b popB -r 1 -o out.tsv \
  -w 200000 -s 100000 -m 2 -t 4
```

## Synopsis

```text
xpclr -i <vcf> -p <pop.txt> -a <popA> -b <popB> -o <out.tsv>
      [-r <region>] [-w <size>] [-s <step>] [-k <maxsnps>] [-m <minsnps>]
      [-L <ld>] [-N <ne>] [-R <rrate>] [-G <gmap>] [--omega-trim <frac>]
      [-t <threads>] [--seed INT] [-P <phased>] [--unimodal-s] [-V INT]
```

### Required

| Option | Description |
|--------|-------------|
| `-i`, `--input` | VCF/BCF (bgzip + TBI/CSI recommended) |
| `-p`, `--pop` | Population file: `SAMPLE  GROUP` |
| `-a`, `--popA` | Target population name |
| `-b`, `--popB` | Reference population name |
| `-o`, `--out` | Output TSV |

### Optional

| Option | Default | Description |
|--------|---------|-------------|
| `-r`, `--regions` | all contigs | Contig or interval (`Chr01`, `Chr01:200-30000`) |
| `-w`, `--size` | 50000 | Window size (bp) |
| `-s`, `--step` | 25000 | Window step (bp) |
| `-k`, `--maxsnps` | 500 | Max SNPs per window (subsample if denser) |
| `-m`, `--minsnps` | 10 | Min SNPs per window (`>= 2`) |
| `-L`, `--ld` | 0.95 | LD r² weight cutoff |
| `-N`, `--ne` | 20000 | Effective population size |
| `-R`, `--rrate` | 1e-8 | Recombination rate per bp without `--gmap` |
| `-G`, `--gmap` | none | Genetic map `CHROM POS GDIST` |
| `--omega-trim` | 0.01 | Drop top fraction of SNP r when estimating ω; `0` = raw mean |
| `-P`, `--phased` | 3 | LD weight mode (see below) |
| `-t`, `--threads` | 1 | Threads |
| `--seed` | 1 | RNG seed for `maxsnps` subsample |
| `--unimodal-s` | off | Stop at first LL decline along s (hardingnj) |
| `-V` | 1 | Log level: 0 quiet, 1 info, 2 debug |

## Input

**VCF/BCF:** diploid `GT`, biallelic SNPs only. Index with `bcftools index file.vcf.gz`.

**Pop file** (`-p`):

```text
# SAMPLE  GROUP
S1  popA
S2  popA
S3  popB
```

Whitespace-separated. `#` comments. Duplicate sample with different groups is an error.

**Genetic map** (`--gmap`, optional): `CHROM POS GDIST`, sorted by POS within chrom. Default genetic distance is `POS * rrate`.

## LD weight mode (`--phased`)

Controls how pairwise SNP correlation (eq 7) is computed from the reference population. Missing genotypes and phasing affect the result.

| Value | Mode | Missing | Notes |
|-------|------|---------|-------------|
| `0` | dosage-fill | fill 0 | Missing treated as ref homozygote. Biased when missing is non-random. |
| `1` | phased haplotype r | skip | Single-haplotype frequencies. Needs phased VCF (`0\|1`). No HWE assumption. |
| `2` | EM two-locus r | skip | EM infers haplotype phase from genotypes. Assumes HWE. |
| `3` | pairwise-complete | skip | Dosage r per SNP pair over shared samples. No HWE assumption. Default. |

## Output

TSV columns (hardingnj-compatible names):

```text
id chrom start stop pos_start pos_stop modelL nullL sel_coef
nSNPs nSNPs_avail xpclr xpclr_norm
```

| Column | Description |
|--------|-------------|
| `modelL` / `nullL` | Best / neutral (s=0) composite log-likelihood |
| `sel_coef` | s on grid at `modelL` |
| `nSNPs` | SNPs used after subsample |
| `nSNPs_avail` | SNPs in window before subsample |
| `xpclr` | `2 * (modelL - nullL)` |
| `xpclr_norm` | z-score of finite `xpclr` values in this run |

ω is estimated once per contig/region. LD weights cost O(k² · n_B) per window (k ≤ `--maxsnps`).

## Test

```bash
make test
```

## Citation

Chen H, Patterson N, Reich D. Population differentiation as a test for selective sweeps. Genome Res. 2010;20(3):393-402. doi:10.1101/gr.100545.109

hardingnj/xpclr: https://github.com/hardingnj/xpclr
