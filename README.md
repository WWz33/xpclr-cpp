# xpclr-cpp

C++ / [htslib](https://github.com/samtools/htslib) rewrite of
[hardingnj/xpclr](https://github.com/hardingnj/xpclr) (Chen, Patterson & Reich 2010).

Goals:

- numerical parity with the Python reference (hardingnj/xpclr)
- faster VCF IO via htslib (+ CSI/TBI region queries)
- OpenMP window-level parallel scan
- population file + GCTA-style sample matching logs
- bioinformatics CLI conventions (`-h`/`-v` exit 0, errors on stderr)

## Build

Dependencies: C++17, OpenMP, htslib, GSL (integration), zlib.

```bash
make -j
./xpclr -h
```

Override htslib if needed:

```bash
make HTS_CFLAGS='-I/path/to/htslib' HTS_LIBS='-L/path/to -lhts'
```

## Usage

```text
xpclr -i <vcf.gz> --pop <pop.txt> -a <popA> -b <popB> -c <chr> -o <out.tsv>
      [--size 20000 --step 20000 --ld 0.95 --rrate 1e-8 --threads N --seed 1]
      [--no-early-stop]
```

### Population file

Two columns (whitespace), optional `#` comments:

```text
# SAMPLE  GROUP
FENGGWS001  popA
FENGGWS002  popA
FENGGWS100  popB
FENGGWS200  popC
```

- `-a` / `-b` select two group names for a single pairwise XP-CLR run
- sample IDs are matched to the VCF header by **name**
- samples not in the selected groups are ignored
- same sample, same group twice: warn, keep once
- same sample, different groups: **error**
- selected group with 0 VCF matches: **error** (exit 1)

### Smoke demo

```bash
# demo/smoke.vcf.gz + demo/pop_smoke.txt are prepared from FENGGWS348_ld0.8.vcf.gz
./xpclr -i demo/smoke.vcf.gz --pop demo/pop_smoke.txt -a popA -b popB -c 1 \
  -o demo/out.tsv --size 200000 --step 100000 --minsnps 2 --threads 4
```

Full 348-sample pop map: `demo/pop.txt` (popA/popB/popC).

## Output

TSV columns aligned with hardingnj/xpclr:

`id chrom start stop pos_start pos_stop modelL nullL sel_coef nSNPs nSNPs_avail xpclr xpclr_norm`

`xpclr = 2 * (modelL - nullL)`, `xpclr_norm` is z-score over finite windows.

## Parity notes

- allele filters: multiallelic / missing-in-pop / popB fixed-or-singleton
- omega, LD weights (Rogers–Huff r² on popB dosage, missing filled as 0)
- selection coefficient grid identical to Python
- default **early-stop** on the selection grid matches Python `compute_xpclr`
- use `--no-early-stop` to evaluate the full grid (see issue #115)
- `--seed` makes `maxsnps` subsampling deterministic (Python uses global numpy RNG)

Smoke vs Python (50 windows, chr1 subset): modelL MAE ~3e-11, sel_coef exact match
with default early-stop.

## Upstream issues (verified)

See [docs/ISSUES.md](docs/ISSUES.md).

| Issue | Verdict | Action in this rewrite |
|-------|---------|------------------------|
| #113 scipy romberg removed | real dependency bug | use GSL QAGS; no scipy |
| #115 / #101 / #107 mostly zero scores | real heuristic (early-stop) | default keeps Python behaviour; `--no-early-stop` available |
| VCF + no `gdistkey` | real (Python builds `variants/None`) | always fall back to `pos * rrate` |
| #105 all multiallelic | input / format, not formula | clear filter counters in log |

## License

Algorithm from Chen et al. 2010; this rewrite is provided for local research use.
Upstream Python package: MIT (hardingnj/xpclr).
