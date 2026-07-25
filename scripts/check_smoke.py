#!/usr/bin/env python3
from pathlib import Path

def rows(path):
    lines = Path(path).read_text().splitlines()
    hdr = lines[0].split("\t")
    out = []
    for line in lines[1:]:
        d = dict(zip(hdr, line.split("\t")))
        out.append(d)
    return out

got = rows("/tmp/xpclr_smoke.tsv")
exp = rows("data/smoke_expected.tsv")
assert len(got) == len(exp), (len(got), len(exp))
keys = ["id", "chrom", "start", "stop", "nSNPs", "nSNPs_avail", "sel_coef"]
for i, (g, e) in enumerate(zip(got, exp)):
    for k in keys:
        assert g[k] == e[k], (i, k, g.get(k), e.get(k))
    # numeric soft check on xpclr when both finite
    if g["xpclr"] != "nan" and e["xpclr"] != "nan":
        gf = float(g["xpclr"]); ef = float(e["xpclr"])
        assert abs(gf - ef) <= max(1e-6, 1e-5 * abs(ef)), (i, gf, ef)
print("test-smoke ok", len(got), "windows")
