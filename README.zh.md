# xpclr-cpp

<!-- README-I18N:START -->

[English](./README.md) | **汉语**

<!-- README-I18N:END -->

基于 **C++ / [htslib](https://github.com/samtools/htslib)** 的 **XP-CLR** 扫描实现（Chen, Patterson & Reich 2010），算法路径对齐常用 Python 包 [hardingnj/xpclr](https://github.com/hardingnj/xpclr)。

## 相对 hardingnj/xpclr 的改动

- **htslib VCF/BCF**：CSI/TBI 区域查询 `chr:start-(stop+size)`；亦可整条 contig 加载
- **群体文件**：一份 `SAMPLE  GROUP` 映射，用 `-a` / `-b` 指定群体
- **默认全 \(s\) 网格最大化**：比 Python early-stop 更少假零（[#115](https://github.com/hardingnj/xpclr/issues/115)）；对齐原版请加 `--unimodal-s`
- **可复现子采样**：`--seed` 控制 `maxsnps` 随机抽薄

## 编译

```bash
git clone <this-repo> xpclr-cpp
cd xpclr-cpp
make -j
./xpclr -h
./xpclr -v
```

若 `pkg-config` 找不到 htslib：

```bash
make HTS_CFLAGS='-I/path/to/include' HTS_LIBS='-L/path/to/lib -lhts'
```

## 快速开始

```bash
bcftools index demo/smoke.vcf.gz
./xpclr -i demo/smoke.vcf.gz -p demo/pop_smoke.txt \
  -a popA -b popB -r 1 -o demo/out.tsv \
  --size 200000 --step 100000 --minsnps 2 --threads 4

bash scripts/prep_smoke.sh /path/to/FENGGWS348_ld0.8.vcf.gz
```

## 用法

```text
xpclr -i <vcf.gz> -p <pop.txt> -a <popA> -b <popB> -o <out.tsv>
      [-r <region>] [--size INT] [--step INT]
      [--maxsnps INT] [--minsnps INT] [--ld FLOAT] [--rrate FLOAT]
      [--threads INT] [--seed INT] [--unimodal-s] [-V INT]
```

### 必选参数

| 参数 | 说明 |
|------|------|
| `-i`, `--input` | VCF/BCF（建议 bgzip + TBI/CSI） |
| `-p`, `--pop` | 群体表：见下方「群体文件」 |
| `-a`, `--popA` | 群体 A 名（选择目标） |
| `-b`, `--popB` | 群体 B 名（参照） |
| `-o`, `--out` | 输出 TSV 路径 |

### 可选参数

| 参数 | 默认 | 说明 |
|------|------|------|
| `-r`, `--regions` | 全部 contig | contig 或区间：`Chr01`、`Chr01:200-30000`、`1:1000000-` |
| `--size` | 20000 | 窗长（bp） |
| `--step` | 20000 | 步长（bp） |
| `--maxsnps` | 200 | 每窗最多 SNP（过密则随机子采样） |
| `--minsnps` | 10 | 每窗最少 SNP（`>= 2`） |
| `--ld` | 0.95 | LD \(r^2\) 权重阈值 |
| `--rrate` | 1e-8 | 无遗传图时的重组率 / bp |
| `--threads` | 1 | 线程数 |
| `--seed` | 1 | `maxsnps` 子采样随机种子 |
| `--unimodal-s` | 关 | 沿 \(s\) 首次似然下降即停（hardingnj/python） |
| `-V` | 1 | 日志：0 quiet，1 info，2 debug |

## 输入

### VCF / BCF

- 二倍体 `GT`；仅保留双等位 SNP（多等位 / indel 丢弃）
- 建议索引：`bcftools index file.vcf.gz`

### 群体文件（`-p` / `--pop`）

空白分隔，索引 Sample Name，无需排序，`#` 开头为注释：

```text
# SAMPLE  GROUP
FENGGWS001  popA
FENGGWS002  popA
FENGGWS100  popB
FENGGWS200  popC
```

## 输出

制表符分隔，列名兼容 hardingnj：

```text
id  chrom  start  stop  pos_start  pos_stop  modelL  nullL  sel_coef
nSNPs  nSNPs_avail  xpclr  xpclr_norm
```

| 列 | 含义 |
|----|------|
| `modelL` / `nullL` | 最优 / 中性（\(s=0\)）复合对数似然 |
| `sel_coef` | 达到 `modelL` 的网格 \(s\) |
| `nSNPs` | 窗内实际使用 SNP（含子采样后） |
| `nSNPs_avail` | 子采样前窗内 SNP 数 |
| `xpclr` | \(2 \times (\mathrm{modelL} - \mathrm{nullL})\) |
| `xpclr_norm` | 本次运行有限窗上的 `xpclr` z 分数 |

## 目录结构

```text
xpclr-cpp/
├── src/           主程序、VCF I/O、群体、XP-CLR 核心
├── include/       xpclr.hpp
├── demo/          smoke VCF、群体表、可选基准输出
├── scripts/       prep_smoke.sh
├── docs/          ISSUES.md
├── Makefile
└── xpclr          make 后生成的可执行文件
```

## 引用

> Chen H, Patterson N, Reich D. Population differentiation as a test for selective sweeps. *Genome Res.* 2010;20(3):393–402. doi:10.1101/gr.100545.109

> [hardingnj/xpclr](https://github.com/hardingnj/xpclr)
