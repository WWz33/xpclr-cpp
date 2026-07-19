# Upstream hardingnj/xpclr issues — verification notes

Checked with `gh search issues --repo hardingnj/xpclr` and source reading of
`methods.py` / `util.py` / `bin/xpclr` (v1.1.2 lineage).

## Confirmed real bugs / design defects

### #113 — scipy>=1.15 removes `romberg`

- **True.** `from scipy.integrate import romberg` fails on SciPy 1.15+.
- Active path uses `quad` (`chen_likelihood`), but import still pulls `romberg`.
- **Fix here:** GSL `gsl_integration_qags` with epsrel=0.001, epsabs=0 (match quad).

### #115 / #101 / #107 — almost all xpclr == 0

- **True as a heuristic defect**, not a segfault.
- `compute_xpclr` walks the selection grid and **breaks on first non-improving ll**:

```python
if ll < maximum_li:
    maximum_li = ll
    maxli_sc = s_coef
else:
    break
```

- Assumes likelihood is unimodal / monotone in the discrete grid. When the first
  step (s=1e-5) is worse than s=0, the scan stops and xpclr stays 0.
- Smoke data: early-stop → 5/50 nonzero windows; full grid → 28/50 nonzero.
- **Fix here:** keep early-stop as default (parity); expose `--no-early-stop`.

### VCF path without genetic map key

- **True.** `load_vcf_wrapper` always requests `variants/{gdistkey}`; default
  `gdistkey=None` becomes field `variants/None` and breaks.
- **Fix here:** no genetic map required; genetic distance = `POS * rrate`.

## Not treated as core formula bugs

| Issue | Notes |
|-------|-------|
| #105 all multiallelic | Usually TXT genotype layout / allele coding; filters work on clean VCF |
| #102 segfault | Environment / scikit-allel / huge windows; not reproduced in C++ path |
| #110 indels/SVs | Method is SNP-oriented; we keep biallelic SNP only |
| usage / interpretation issues | Docs / biology, not code defects |

## Recommendation

For production scans on diverged pops where many zeros look wrong, run with
`--no-early-stop` and compare. Report both if publishing methods.
