#!/usr/bin/env bash
# Build smoke VCF + 3-pop file from FENGGWS348_ld0.8.vcf.gz
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
VCF_SRC="${1:-/home/ww/xpclr/FENGGWS348_ld0.8.vcf.gz}"
OUTDIR="$ROOT/demo"
mkdir -p "$OUTDIR"

echo "[prep] indexing source if needed..."
if [[ ! -f "${VCF_SRC}.tbi" && ! -f "${VCF_SRC}.csi" ]]; then
  bcftools index -t -f "$VCF_SRC"
fi

echo "[prep] full 3-pop map..."
mapfile -t S < <(bcftools query -l "$VCF_SRC")
N=${#S[@]}
A=$((N/3)); B=$((2*N/3))
{
  echo "# SAMPLE GROUP  (auto from $VCF_SRC)"
  for ((i=0;i<A;i++)); do echo -e "${S[i]}\tpopA"; done
  for ((i=A;i<B;i++)); do echo -e "${S[i]}\tpopB"; done
  for ((i=B;i<N;i++)); do echo -e "${S[i]}\tpopC"; done
} > "$OUTDIR/pop.txt"

{
  echo "# smoke 3 pops x 15"
  for ((i=0;i<15;i++)); do echo -e "${S[i]}\tpopA"; done
  for ((i=15;i<30;i++)); do echo -e "${S[i]}\tpopB"; done
  for ((i=30;i<45;i++)); do echo -e "${S[i]}\tpopC"; done
} > "$OUTDIR/pop_smoke.txt"

echo "[prep] smoke VCF chr1:1-5Mb, 45 samples..."
KEEP=$(printf '%s,' "${S[@]:0:45}"); KEEP=${KEEP%,}
bcftools view -r 1:1-5000000 -Ou "$VCF_SRC" \
  | bcftools view -s "$KEEP" -Oz -o "$OUTDIR/smoke.vcf.gz"
bcftools index -t -f "$OUTDIR/smoke.vcf.gz"
echo "[prep] sites: $(bcftools view -H "$OUTDIR/smoke.vcf.gz" | wc -l)"
echo "done: $OUTDIR/smoke.vcf.gz $OUTDIR/pop.txt $OUTDIR/pop_smoke.txt"
