#!/usr/bin/env python3
import re
from pathlib import Path

def omega(path):
    text = Path(path).read_text()
    m = re.search(r"Omega estimated as : ([0-9.eE+-]+)", text)
    assert m, path
    return float(m.group(1)), text

tr, tlog = omega("/tmp/xpclr_omega_trim.log")
rw, rlog = omega("/tmp/xpclr_omega_raw.log")
assert "trim=0.010000" in tlog
assert "trim=0.000000" in rlog
assert rw > 0 and tr <= rw + 1e-12, (tr, rw)
print("test-omega ok", "raw", rw, "trim", tr)
