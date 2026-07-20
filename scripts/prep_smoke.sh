#!/usr/bin/env bash
# Build data/smoke.vcf.gz + data/pop_smoke.txt from a large VCF
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
VCF_SRC="${1:-/home/ww/xpclr/FENGGWS348_ld0.8.vcf.gz}"
OUTDIR="$ROOT/data"
mkdir -p "$OUTDIR"

if [[ ! -f "${VCF_SRC}.tbi" && ! -f "${VCF_SRC}.csi" ]]; then
  bcftools index -t -f "$VCF_SRC"
fi

mapfile -t S < <(bcftools query -l "$VCF_SRC")
{
  echo "# SAMPLE GROUP"
  for ((i=0;i<15;i++)); do echo -e "${S[i]}\tpopA"; done
  for ((i=15;i<30;i++)); do echo -e "${S[i]}\tpopB"; done
  for ((i=30;i<45;i++)); do echo -e "${S[i]}\tpopC"; done
} > "$OUTDIR/pop_smoke.txt"

KEEP=$(printf '%s,' "${S[@]:0:45}"); KEEP=${KEEP%,}
bcftools view -r 1:1-5000000 -Ou "$VCF_SRC" \
  | bcftools view -s "$KEEP" -Oz -o "$OUTDIR/smoke.vcf.gz"
bcftools index -t -f "$OUTDIR/smoke.vcf.gz"
echo "done: $OUTDIR/smoke.vcf.gz $OUTDIR/pop_smoke.txt"
