# xpclr-cpp

<!-- README-I18N:START -->

[English](./README.md) | **汉语**

<!-- README-I18N:END -->

基于 **C++ / [htslib](https://github.com/samtools/htslib)** 的 **XP-CLR** 扫描实现（Chen, Patterson & Reich 2010），算法路径对齐常用 Python 包 [hardingnj/xpclr](https://github.com/hardingnj/xpclr)。

特点：更快的 VCF I/O（索引区域拉取 + BGZF 多线程）、OpenMP 窗并行、单文件群体表、默认在选择系数 \(s\) 上 **全网格取最大**（可用 `--unimodal-s` 做 hardingnj 风格提前退出）。

方法原文：[Genome Res. 2010](https://www.ncbi.nlm.nih.gov/pubmed/20086244) · 上游 issue 核对：[docs/ISSUES.md](docs/ISSUES.md)

## 相对 hardingnj/xpclr 的改动

- **htslib VCF/BCF**：CSI/TBI 区域查询 `chr:start-(stop+size)`；亦可整条 contig 加载
- **OpenMP**：`--threads N` 同时用于 BGZF 解压与窗扫描
- **群体文件**：一份 `SAMPLE  GROUP` 映射，用 `-a` / `-b` 点名两群（GCTA 风格 input/matched 日志）
- **默认全 \(s\) 网格最大化**：比 Python early-stop 更少假零（[#115](https://github.com/hardingnj/xpclr/issues/115)）；对齐原版请加 `--unimodal-s`
- **可复现子采样**：`--seed` 控制 `maxsnps` 随机抽薄
- **CLI 惯例**：`-h` / `-v` 退出码 0；错误走 stderr 并提示 `-h`

## 依赖

| 组件 | 用途 |
|------|------|
| C++17 + OpenMP | 编译 / 窗并行 |
| [htslib](https://github.com/samtools/htslib) | VCF/BCF I/O |
| [GSL](https://www.gnu.org/software/gsl/) | 自适应积分（对应 scipy.quad） |
| zlib, libbz2, liblzma, libcurl | htslib 链接依赖 |

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
# smoke：小片段 chr1 + 三个群体
./xpclr -i demo/smoke.vcf.gz --pop demo/pop_smoke.txt \
  -a popA -b popB -c 1 -o demo/out.tsv \
  --size 200000 --step 100000 --minsnps 2 --threads 4

# 从大 VCF 重建 smoke（可选）
bash scripts/prep_smoke.sh /path/to/FENGGWS348_ld0.8.vcf.gz
```

## 用法

```text
xpclr -i <vcf.gz> --pop <pop.txt> -a <popA> -b <popB> -c <chr> -o <out.tsv>
      [--size INT] [--step INT] [--start INT] [--stop INT]
      [--maxsnps INT] [--minsnps INT] [--ld FLOAT] [--rrate FLOAT]
      [--threads INT] [--seed INT] [--unimodal-s] [-V INT]
```

### 必选参数

| 参数 | 说明 |
|------|------|
| `-i`, `--input` | VCF/BCF（建议 bgzip + TBI/CSI） |
| `--pop` | 群体表：两列 `SAMPLE  GROUP` |
| `-a`, `--popA` | 群体 A 名（选择目标） |
| `-b`, `--popB` | 群体 B 名（参照） |
| `-c`, `--chr` | contig 名，**须与 VCF header 完全一致**（如 `1` 或 `Chr01`） |
| `-o`, `--out` | 输出 TSV 路径 |

### 常用选项

| 参数 | 默认 | 说明 |
|------|------|------|
| `--size` | 20000 | 窗长（bp） |
| `--step` | 20000 | 步长（bp） |
| `--start` | 1 | 首窗起点；同时限制 VCF 加载下界 |
| `--stop` | 0 | 窗扫描上界；`0` = 加载区间末 SNP（不写 stop 则整 contig） |
| `--maxsnps` | 200 | 每窗最多 SNP（过密则随机子采样） |
| `--minsnps` | 10 | 每窗最少 SNP（`>= 2`） |
| `--ld` | 0.95 | LD \(r^2\) 权重阈值 |
| `--rrate` | 1e-8 | 无遗传图时的重组率 / bp |
| `--threads` | 1 | htslib BGZF 线程 + OpenMP 窗并行 |
| `--seed` | 1 | `maxsnps` 子采样随机种子 |
| `--unimodal-s` | 关 | 沿 \(s\) 首次似然下降即停（hardingnj/python） |
| `-V` | 1 | 日志：0 quiet，1 info，2 debug |

PBS 风格示例（密窗、10 线程）：

```bash
./xpclr -i FENGGWS348.vcf.gz --pop pops.txt -a W -b C -c Chr01 \
  -o W_vs_C.Chr01.xpclr.tsv \
  --size 20000 --step 2000 --maxsnps 300 --minsnps 10 \
  --threads 10 --seed 1
```

## 输入

### VCF / BCF

- 二倍体 `GT`；仅保留双等位 SNP（多等位 / indel 丢弃）
- 建议索引：`bcftools index -t file.vcf.gz`
- `-c` 必须与 header contig 字符串一致

### 群体文件（`--pop`）

空白分隔，`#` 开头为注释：

```text
# SAMPLE  GROUP
FENGGWS001  popA
FENGGWS002  popA
FENGGWS100  popB
FENGGWS200  popC
```

| 规则 | 行为 |
|------|------|
| Sample ID | 与 VCF header sample 名匹配 |
| `-a` / `-b` | 本次只用这两个 group；其余 group 忽略 |
| 同 sample 同 group 重复 | 警告，保留一次 |
| 同 sample 不同 group | **报错**，exit 1 |
| 选中 group 在 VCF 中 0 匹配 | **报错**，exit 1 |
| 在 VCF 但不在 `-a`/`-b` | 忽略 |

日志打印每个选中 group 的 **input N / matched / used**（GCTA 风格）。

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

## 说明

### I/O 与区域

- `--threads N`：加载阶段 htslib BGZF 线程，随后 OpenMP 扫窗
- `--stop > 0` 时加载 `chr:start-(stop+size)`，保证最后一窗仍有 SNP
- **Omega** 仅在**已加载** SNP 上估计。区域加载与整 chrom 加载的分数可不同
- 与 hardingnj 数值对齐：整 chrom（`--stop 0`）+ `--unimodal-s`

### 过滤（与 hardingnj 一致思路）

1. 多等位 / 非 SNP  
2. 任一族群全 missing  
3. 群体 B 固定或 singleton  

### 选择系数网格

默认：在完整 \(s\) 网格  
`{0, 1e-5, 5e-5, …, 0.15}` 上最大化复合似然（选择扫描推荐）。

`--unimodal-s`：沿网格首次似然下降即停（Python `compute_xpclr`；可能抬高零分窗，见 [docs/ISSUES.md](docs/ISSUES.md)）。

### 与 Python 对齐检查表

| 需求 | 设置 |
|------|------|
| 与 Python 同 chrom SNP 集合 | `--stop 0` |
| 与 Python 同 early-exit | `--unimodal-s` |
| maxsnps 子采样 | Python 无固定 seed；C++ 用 `--seed`（触发子采样时期望不完全一致） |

Smoke（50 窗，`demo/smoke.vcf.gz`）+ `--unimodal-s`：modelL MAE \(\sim 10^{-11}\)，`sel_coef` 完全一致。

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

使用 XP-CLR 方法请引用：

> Chen H, Patterson N, Reich D. Population differentiation as a test for selective sweeps. *Genome Res.* 2010;20(3):393–402. doi:10.1101/gr.100545.109

与 hardingnj 结果对比时请同时引用：

> [hardingnj/xpclr](https://github.com/hardingnj/xpclr)

## 许可

算法：Chen et al. 2010。上游 hardingnj/xpclr 为 MIT。本重写供研究使用；公开发布前请补充项目 `LICENSE`。
