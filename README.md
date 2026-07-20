# xpclr-cpp

<!-- README-I18N:START -->

**English** | [汉语](./README.zh.md)

<!-- README-I18N:END -->

C++ / [htslib](https://github.com/samtools/htslib) implementation of the **XP-CLR** scan (Chen, Patterson & Reich 2010), aligned with the widely used Python package [hardingnj/xpclr](https://github.com/hardingnj/xpclr).

Faster VCF I/O (indexed region fetch + BGZF threads), OpenMP window parallelisation, a single population map file, and a default **full-grid** maximisation over the selection coefficient \(s\) (optional hardingnj-style early exit via `--unimodal-s`).

Method paper: [Genome Res. 2010](https://www.ncbi.nlm.nih.gov/pubmed/20086244) · Upstream issues notes: [docs/ISSUES.md](docs/ISSUES.md)

## Features (vs hardingnj/xpclr)

- **htslib VCF/BCF** — CSI/TBI region query `chr:start-(stop+size)`; optional whole-contig load
- **OpenMP** — `--threads N` for BGZF decompression and window scan
- **Population file** — one `SAMPLE  GROUP` map; select two groups with `-a` / `-b` (GCTA-style match logs)
- **Default full \(s\)-grid max** — fewer false-zero windows than Python early-stop ([#115](https://github.com/hardingnj/xpclr/issues/115)); use `--unimodal-s` for hardingnj parity
- **Reproducible subsample** — `--seed` for `maxsnps` random thinning
- **CLI conventions** — `-h` / `-v` exit 0; errors on stderr; points to `-h`

## Dependencies

| Component | Role |
|-----------|------|
| C++17 compiler + OpenMP | build / window parallel |
| [htslib](https://github.com/samtools/htslib) | VCF/BCF I/O |
| [GSL](https://www.gnu.org/software/gsl/) | adaptive integration (scipy.quad analogue) |
| zlib, libbz2, liblzma, libcurl | htslib link deps |

## Build

```bash
git clone <this-repo> xpclr-cpp
cd xpclr-cpp
make -j
./xpclr -h
./xpclr -v
```

If htslib is not found via `pkg-config`:

```bash
make HTS_CFLAGS='-I/path/to/include' HTS_LIBS='-L/path/to/lib -lhts'
```

## Quick start

```bash
# smoke demo (small chr1 subset + 3 groups)
./xpclr -i demo/smoke.vcf.gz --pop demo/pop_smoke.txt \
  -a popA -b popB -c 1 -o demo/out.tsv \
  --size 200000 --step 100000 --minsnps 2 --threads 4

# rebuild smoke inputs from a large VCF (optional)
bash scripts/prep_smoke.sh /path/to/FENGGWS348_ld0.8.vcf.gz
```

## Usage

```text
xpclr -i <vcf.gz> --pop <pop.txt> -a <popA> -b <popB> -c <chr> -o <out.tsv>
      [--size INT] [--step INT] [--start INT] [--stop INT]
      [--maxsnps INT] [--minsnps INT] [--ld FLOAT] [--rrate FLOAT]
      [--threads INT] [--seed INT] [--unimodal-s] [-V INT]
```

### Required

| Flag | Description |
|------|-------------|
| `-i`, `--input` | VCF/BCF (bgzip + TBI/CSI recommended) |
| `--pop` | Population map: two columns `SAMPLE  GROUP` |
| `-a`, `--popA` | Group name for population A (selected / target) |
| `-b`, `--popB` | Group name for population B (reference) |
| `-c`, `--chr` | Contig name **exactly as in the VCF header** (e.g. `1` or `Chr01`) |
| `-o`, `--out` | Output TSV path |

### Common options

| Flag | Default | Description |
|------|---------|-------------|
| `--size` | 20000 | Window size (bp) |
| `--step` | 20000 | Window step (bp) |
| `--start` | 1 | First window start; also trims VCF load lower bound |
| `--stop` | 0 | Last base for windows; `0` = last loaded SNP (whole contig if no stop) |
| `--maxsnps` | 200 | Max SNPs per window (random subsample if denser) |
| `--minsnps` | 10 | Min SNPs per window (`>= 2`) |
| `--ld` | 0.95 | LD \(r^2\) cutoff for SNP weights |
| `--rrate` | 1e-8 | Recombination rate per bp if no genetic map |
| `--threads` | 1 | htslib BGZF threads + OpenMP windows |
| `--seed` | 1 | RNG seed for `maxsnps` subsample |
| `--unimodal-s` | off | Stop at first likelihood decline along \(s\) (hardingnj/python) |
| `-V` | 1 | Log level: 0 quiet, 1 info, 2 debug |

PBS-style example (dense windows, 10 threads):

```bash
./xpclr -i FENGGWS348.vcf.gz --pop pops.txt -a W -b C -c Chr01 \
  -o W_vs_C.Chr01.xpclr.tsv \
  --size 20000 --step 2000 --maxsnps 300 --minsnps 10 \
  --threads 10 --seed 1
```

## Input

### VCF / BCF

- Diploid `GT` field; biallelic SNPs only (multiallelic / indels dropped).
- Index recommended: `bcftools index -t file.vcf.gz`
- Contig id must match `-c` exactly.

### Population file (`--pop`)

Whitespace-separated, optional `#` comments:

```text
# SAMPLE  GROUP
FENGGWS001  popA
FENGGWS002  popA
FENGGWS100  popB
FENGGWS200  popC
```

| Rule | Behaviour |
|------|-----------|
| Sample ID | Match VCF header sample names |
| `-a` / `-b` | Two group names used in this run; other groups ignored |
| Same sample, same group twice | warn, keep once |
| Same sample, different groups | **error**, exit 1 |
| Selected group with 0 VCF matches | **error**, exit 1 |
| Sample in VCF but not in `-a`/`-b` | ignored |

Log lines report **input N / matched / used** per selected group (GCTA-style).

## Output

Tab-separated columns (hardingnj-compatible names):

```text
id  chrom  start  stop  pos_start  pos_stop  modelL  nullL  sel_coef
nSNPs  nSNPs_avail  xpclr  xpclr_norm
```

| Column | Meaning |
|--------|---------|
| `modelL` / `nullL` | Best / neutral (\(s=0\)) composite log-likelihood |
| `sel_coef` | Grid \(s\) attaining `modelL` |
| `nSNPs` | SNPs used in window (after optional subsample) |
| `nSNPs_avail` | SNPs in window before subsample |
| `xpclr` | \(2 \times (\mathrm{modelL} - \mathrm{nullL})\) |
| `xpclr_norm` | Z-score of `xpclr` over finite windows on this run |

## Notes

### I/O and regions

- `--threads N` sets htslib BGZF threads during load, then OpenMP over windows.
- With `--stop > 0`, variants are loaded as `chr:start-(stop+size)` so the last window still sees SNPs.
- **Omega** (cross-population variance scalar) is estimated on **loaded** SNPs only. Regional load can change scores vs whole-contig load even on shared windows.
- For hardingnj numerical parity: whole-contig load (`--stop 0`) + `--unimodal-s`.

### Filters (same spirit as hardingnj)

1. Multiallelic / non-SNP  
2. All missing in either population  
3. Invariant or singleton in population B  

### Selection grid

Default: maximise composite likelihood over the full \(s\) grid  
`{0, 1e-5, 5e-5, …, 0.15}` (recommended for selection scans).

`--unimodal-s`: stop at the first likelihood decline along that grid (Python `compute_xpclr` behaviour; can inflate zero scores — see [docs/ISSUES.md](docs/ISSUES.md)).

### Parity checklist vs Python

| Requirement | Flag / setting |
|-------------|----------------|
| Same SNP set as Python chrom load | `--stop 0` (whole contig) |
| Same early-exit heuristic | `--unimodal-s` |
| Same maxsnps draws | Python has no fixed seed; C++ uses `--seed` (expect mismatch when subsample fires) |

Smoke (50 windows, `demo/smoke.vcf.gz`) with `--unimodal-s`: modelL MAE \(\sim 10^{-11}\), `sel_coef` exact match.

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

If you use XP-CLR, cite the method:

> Chen H, Patterson N, Reich D. Population differentiation as a test for selective sweeps. *Genome Res.* 2010;20(3):393–402. doi:10.1101/gr.100545.109

Also cite the Python reference implementation when comparing results:

> [hardingnj/xpclr](https://github.com/hardingnj/xpclr)

## License

Algorithm: Chen et al. 2010. Upstream hardingnj/xpclr is MIT. This rewrite is provided for research use; add a project `LICENSE` before public redistribution.
