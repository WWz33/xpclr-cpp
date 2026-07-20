# xpclr-cpp

<!-- README-I18N:START -->

**English** | [汉语](./README.zh.md)

<!-- README-I18N:END -->

C++ / [htslib](https://github.com/samtools/htslib) XP-CLR implementation (Chen, Patterson & Reich 2010), algorithm path aligned with [hardingnj/xpclr](https://github.com/hardingnj/xpclr).

## Changes vs hardingnj/xpclr

- **htslib VCF/BCF** — CSI/TBI region query `chr:start-(stop+size)`; whole-contig load also supported
- **Population file** — one `SAMPLE  GROUP` map; select groups with `-a` / `-b`
- **Default full \(s\)-grid max** — fewer false zeros than Python early-stop ([#115](https://github.com/hardingnj/xpclr/issues/115)); use `--unimodal-s` for hardingnj parity
- **Reproducible subsample** — `--seed` for `maxsnps` thinning

## Build

```bash
git clone <this-repo> xpclr-cpp
cd xpclr-cpp
make -j
./xpclr -h
./xpclr -v
```

If `pkg-config` cannot find htslib:

```bash
make HTS_CFLAGS='-I/path/to/include' HTS_LIBS='-L/path/to/lib -lhts'
```

## Quick start

```bash
bcftools index demo/smoke.vcf.gz
./xpclr -i demo/smoke.vcf.gz -p demo/pop_smoke.txt \
  -a popA -b popB -r 1 -o demo/out.tsv \
  --size 200000 --step 100000 --minsnps 2 --threads 4

bash scripts/prep_smoke.sh /path/to/FENGGWS348_ld0.8.vcf.gz
```

## Usage

```text
xpclr -i <vcf.gz> -p <pop.txt> -a <popA> -b <popB> -o <out.tsv>
      [-r <region>] [--size INT] [--step INT]
      [--maxsnps INT] [--minsnps INT] [--ld FLOAT] [--rrate FLOAT]
      [--threads INT] [--seed INT] [--unimodal-s] [-V INT]
```

### Required

| Flag | Description |
|------|-------------|
| `-i`, `--input` | VCF/BCF (bgzip + TBI/CSI recommended) |
| `-p`, `--pop` | Population map (see **Population file** below) |
| `-a`, `--popA` | Population A name (selection target) |
| `-b`, `--popB` | Population B name (reference) |
| `-o`, `--out` | Output TSV path |

### Optional

| Flag | Default | Description |
|------|---------|-------------|
| `-r`, `--regions` | all contigs | Contig or interval: `Chr01`, `Chr01:200-30000`, `1:1000000-` |
| `--size` | 20000 | Window size (bp) |
| `--step` | 20000 | Window step (bp) |
| `--maxsnps` | 200 | Max SNPs per window (random subsample if denser) |
| `--minsnps` | 10 | Min SNPs per window (`>= 2`) |
| `--ld` | 0.95 | LD \(r^2\) weight cutoff |
| `--rrate` | 1e-8 | Recombination rate per bp if no genetic map |
| `--threads` | 1 | Number of threads |
| `--seed` | 1 | RNG seed for `maxsnps` subsample |
| `--unimodal-s` | off | Stop at first likelihood decline along \(s\) (hardingnj/python) |
| `-V` | 1 | Log level: 0 quiet, 1 info, 2 debug |

## Input

### VCF / BCF

- Diploid `GT`; biallelic SNPs only (multiallelic / indels dropped)
- Index recommended: `bcftools index file.vcf.gz`

### Population file (`-p` / `--pop`)

Whitespace-separated; keyed by sample name; order free; `#` starts a comment:

```text
# SAMPLE  GROUP
FENGGWS001  popA
FENGGWS002  popA
FENGGWS100  popB
FENGGWS200  popC
```

## Output

Tab-separated; column names compatible with hardingnj:

```text
id  chrom  start  stop  pos_start  pos_stop  modelL  nullL  sel_coef
nSNPs  nSNPs_avail  xpclr  xpclr_norm
```

| Column | Meaning |
|--------|---------|
| `modelL` / `nullL` | Best / neutral (\(s=0\)) composite log-likelihood |
| `sel_coef` | Grid \(s\) attaining `modelL` |
| `nSNPs` | SNPs used in window (after subsample if any) |
| `nSNPs_avail` | SNPs in window before subsample |
| `xpclr` | \(2 \times (\mathrm{modelL} - \mathrm{nullL})\) |
| `xpclr_norm` | Z-score of `xpclr` over finite windows in this run |

## Layout

```text
xpclr-cpp/
├── src/           main, VCF I/O, pop, XP-CLR core
├── include/       xpclr.hpp
├── demo/          smoke VCF, pop maps, optional benches
├── scripts/       prep_smoke.sh
├── docs/          ISSUES.md
├── Makefile
└── xpclr          binary after make
```

## Citation

> Chen H, Patterson N, Reich D. Population differentiation as a test for selective sweeps. *Genome Res.* 2010;20(3):393–402. doi:10.1101/gr.100545.109

> [hardingnj/xpclr](https://github.com/hardingnj/xpclr)
