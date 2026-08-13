# xpclr-cpp

<!-- README-I18N:START -->

[English](./README.md) | **汉语**

<!-- README-I18N:END -->

C++ XP-CLR（Chen, Patterson & Reich 2010），htslib 读写 VCF/BCF。

相对 [hardingnj/xpclr](https://github.com/hardingnj/xpclr)：

- TBI/CSI 区域查询（`-r`）
- 单文件 `SAMPLE GROUP` 群体表（`-a` / `-b`）
- 默认全 s 网格（`--unimodal-s` 对齐 Python early-stop）
- `--seed` 控制 `maxsnps` 子采样；`--omega-trim` 估计背景 ω；可选 `--gmap`
- `--phased INT` 选择 LD 权重的缺失处理与相位推断方式

## 安装

```bash
git clone --recurse-submodules https://github.com/WWz33/xpclr-cpp.git
cd xpclr-cpp
make -j
```

需要 C++17、zlib、bzip2、lzma、libcurl、openssl、libdeflate。
默认编译 vendored htslib 与 GSL。系统库：`make USE_SYSTEM_HTS=1 USE_SYSTEM_GSL=1`。

## 示例

```bash
./xpclr -i data/smoke.vcf.gz -p data/pop_smoke.txt \
  -a popA -b popB -r 1 -o out.tsv \
  --size 200000 --step 100000 --minsnps 2 --threads 4
```

## 用法

```text
xpclr -i <vcf> -p <pop.txt> -a <popA> -b <popB> -o <out.tsv>
      [-r <region>] [--size INT] [--step INT] [--maxsnps INT] [--minsnps INT]
      [--ld FLOAT] [--rrate FLOAT] [--gmap FILE] [--omega-trim FLOAT]
      [--threads INT] [--seed INT] [--phased INT] [--unimodal-s] [-V INT]
```

### 必选

| 参数 | 说明 |
|------|------|
| `-i`, `--input` | VCF/BCF（建议 bgzip + TBI/CSI） |
| `-p`, `--pop` | 群体文件：`SAMPLE  GROUP` |
| `-a`, `--popA` | 目标群体名 |
| `-b`, `--popB` | 参照群体名 |
| `-o`, `--out` | 输出 TSV |

### 可选

| 参数 | 默认 | 说明 |
|------|------|------|
| `-r`, `--regions` | 全部 contig | contig 或区间（`Chr01`、`Chr01:200-30000`） |
| `--size` | 20000 | 窗长（bp） |
| `--step` | 20000 | 步长（bp） |
| `--maxsnps` | 200 | 每窗最多 SNP（过密则子采样） |
| `--minsnps` | 10 | 每窗最少 SNP（`>= 2`） |
| `--ld` | 0.95 | LD r² 权重阈值 |
| `--rrate` | 1e-8 | 无 `--gmap` 时重组率 / bp |
| `--gmap` | 无 | 遗传图 `CHROM POS GDIST` |
| `--omega-trim` | 0.01 | 估计 ω 时丢掉最高比例的 SNP r；`0` = 原始均值 |
| `--phased` | 3 | LD 权重模式（见下） |
| `--threads` | 1 | 线程数 |
| `--seed` | 1 | `maxsnps` 子采样种子 |
| `--unimodal-s` | 关 | 沿 s 首次似然下降即停（hardingnj） |
| `-V` | 1 | 日志：0 quiet，1 info，2 debug |

## 输入

**VCF/BCF：** 二倍体 `GT`，仅双等位 SNP。索引：`bcftools index file.vcf.gz`。

**群体文件**（`-p`）：

```text
# SAMPLE  GROUP
S1  popA
S2  popA
S3  popB
```

空白分隔。`#` 为注释。同一样本分到不同群体会报错。

**遗传图**（`--gmap`，可选）：`CHROM POS GDIST`，同 chrom 内 POS 升序。默认遗传距离为 `POS * rrate`。

## LD 权重模式（`--phased`）

决定参照群体的成对 SNP 相关（论文 eq 7）怎么算。缺失基因型和相位会影响结果。

| 值 | 模式 | 缺失 | 说明 |
|----|------|------|------------|
| `0` | dosage-fill | 填 0 | 缺失当参考纯合。缺失非随机时偏差大。 |
| `1` | phased 单倍型 r | 跳过 | 直接用单倍型频率算 r。需 phased VCF（`0\|1`）。不假设 HWE。 |
| `2` | EM 两基因型相位 r | 跳过 | EM 从基因型推断单倍型相位。假设 HWE。 |
| `3` | pairwise-complete | 跳过 | 每对 SNP 只用两者都有基因型的样本算 dosage r。不假设 HWE。默认。 |

## 输出

TSV 列名兼容 hardingnj：

```text
id chrom start stop pos_start pos_stop modelL nullL sel_coef
nSNPs nSNPs_avail xpclr xpclr_norm
```

| 列 | 说明 |
|----|------|
| `modelL` / `nullL` | 最优 / 中性（s=0）复合对数似然 |
| `sel_coef` | 达到 `modelL` 的网格 s |
| `nSNPs` | 子采样后使用的 SNP 数 |
| `nSNPs_avail` | 子采样前窗内 SNP 数 |
| `xpclr` | `2 * (modelL - nullL)` |
| `xpclr_norm` | 本次运行有限 `xpclr` 的 z 分数 |

ω 每个 contig/region 估计一次。每窗 LD 权重 O(k² · n_B)（k ≤ `--maxsnps`）。

## 测试

```bash
make test
```

## 引用

Chen H, Patterson N, Reich D. Population differentiation as a test for selective sweeps. Genome Res. 2010;20(3):393-402. doi:10.1101/gr.100545.109

hardingnj/xpclr: https://github.com/hardingnj/xpclr
