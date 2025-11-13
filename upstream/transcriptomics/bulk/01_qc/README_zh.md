# 原始数据质控（Quality Control, QC）
在 RNA-seq 全流程分析中，原始测序数据的质量控制（QC） 是最重要的上游步骤之一。它决定了下游对齐、定量和差异分析的可靠性。原始测序文件通常以 FASTQ 格式 存储，包含每条测序序列（read）及其对应的碱基质量评分。通过对这些原始文件进行 QC，可以发现潜在的问题（如低质量碱基、接头污染、GC 偏倚等），从而决定是否进行修剪或重新测序。

## 1.什么是原始数据质控？
FastQC 是最常用的高通量测序质量评估工具之一，能够快速生成 HTML 报告，直观显示测序质量的多种指标。
### FastQC的主要检测项目包括：

| 模块                           | 说明              | 常见异常原因       |
| ---------------------------- | --------------- | ------------ |
| Per base sequence quality    | 各碱基位置的质量分布      | 测序末端退化、测序仪问题 |
| Per sequence quality scores  | 各 read 的平均质量分布  | 个别低质量 reads  |
| Per base sequence content    | 各位置 A/T/G/C 比例  | 文库建库偏倚、引物污染  |
| Per sequence GC content      | 各 reads 的 GC 含量 | 样本偏倚或污染      |
| Per base N content           | 含 N 的比例         | 测序读段不确定性     |
| Sequence length distribution | read 长度分布       | 过滤/修剪问题      |
| Adapter content              | 接头污染检测          | 文库连接过度、未修剪   |
| Overrepresented sequences    | 过度重复序列          | rRNA 污染、接头序列 |
| K-mer content                | 重复 k-mer 模式     | 特定序列污染或引物残留  |

> FastQC 会为每个模块打分（Pass / Warn / Fail），可快速定位问题。

---

## 2.常用工具与依赖
上游的处理一般都直接在Linux系统中用命令行形式处理，所以常用工具与依赖大部分均在Linux系统中。

### 常用工具
- **FastQC**：质控报告生成（核心工具）
- **MultiQC**：聚合多个样本的 FastQC 报告，便于批量比较
- **TrimGalore / fastp / Cutadapt**：根据 QC 结果进行接头去除与低质量过滤（QC 后处理）

### 安装方法
1.在conda中安装
```bash
conda install -c bioconda fastqc multiqc
```

2.直接下载可执行文件
```bash
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
unzip fastqc_v0.12.1.zip
chmod +x FastQC/fastqc
```

## 3.数据输入与准备
### 输入数据
输入的数据为`.fastq`或`.fastq.gz`，内容为：
```
@序列ID
ACGT...
+
IIIII...
```
含义：

### 数据组织建议
由于上游分析分析步骤较多，每一步会产生很多的文件，所以需要管理好整个项目，推荐一下的文件树管理：
```bash
project/
│── raw_data/          # 原始FASTQ文件
│    ├── sample1_R1.fastq.gz
│    ├── sample1_R2.fastq.gz
│── qc_results/        # 质控输出
│── trimmed_data/      # 修剪后FASTQ（后续步骤）
```

---

## 4.分析流程
### 4.1 单样本质控（FastQC）
```bash
fastqc raw_data/sample1_R1.fastq.gz -o qc_results/
fastqc raw_data/sample1_R2.fastq.gz -o qc_results/
```
执行后将在 `qc_results/` 目录中生成：
- `.html` 报告（可直接浏览器查看）
- `.zip` 结果文件（包含详细统计数据）

### 4.2 批量质控
```bash
mkdir -p qc_results
fastqc raw_data/*.fastq.gz -o qc_results/ -t 4
```
注意：`-t` 参数指定并行线程数，加速批量分析。

### 4.3 聚合报告（MultiQC）
FastQC 为每个文件生成单独报告，不便于批量比较。使用 MultiQC 可以整合多个 FastQC 结果为一个交互式 HTML 报告：
```bash
multiqc qc_results/ -o qc_results/
```
输出：
- `multiqc_report.html`（可交互查看所有样本质量趋势）
- `multiqc_data/`文件夹（包含解析后的统计文件）

---

## 5.结果解释
### 5.1 FastQC报告结构
典型 FastQC HTML 报告由若干模块组成，每个模块会标记为：
- ✅ PASS：质量良好
- ⚠️ WARN：存在轻微异常
- ❌ FAIL：需要关注

### 5.2 关键图表解读
| 模块                            | 理想状态         | 异常信号        | 处理建议                    |
| ----------------------------- | ------------ | ----------- | ----------------------- |
| **Per Base Quality**          | Q30 以上（绿色区域） | 后端下降（黄色/红色） | 修剪低质量末端                 |
| **GC Content**                | 与物种参考分布相似    | 双峰或极端偏移     | 检查污染或物种错误               |
| **Adapter Content**           | 接头含量接近0      | 接头峰上升       | 执行去接头（TrimGalore/fastp） |
| **Sequence Length**           | 一致分布         | 长度不一        | 检查是否混合文库                |
| **Overrepresented Sequences** | None         | 出现>1%重复     | 检查 rRNA/adapter 序列      |
> **一般来说，FASTQ 的整体质量曲线若大部分碱基在 Q30 以上且接头污染少，即可认为质量合格**。


## 6.下游处理（根据 QC 结果）
### 6.1 接头与低质量修剪
若检测到明显的低质量尾端或接头污染，需使用 TrimGalore 或 fastp 进行修剪。参考修剪命令如下：
```bash
# TrimGalore（基于 cutadapt）
trim_galore --paired raw_data/sample1_R1.fastq.gz raw_data/sample1_R2.fastq.gz -o trimmed_data/

# 或使用 fastp（更快）
fastp -i raw_data/sample1_R1.fastq.gz -I raw_data/sample1_R2.fastq.gz \
      -o trimmed_data/sample1_R1.trimmed.fastq.gz \
      -O trimmed_data/sample1_R2.trimmed.fastq.gz \
      --detect_adapter_for_pe --thread 4 --html trimmed_data/sample1_fastp.html
```
**注意：修剪后再运行一次 FastQC 检查修剪效果**。

### 6.2 样本对比
通过 MultiQC 报告可对比不同样本的质量趋势，识别潜在异常样本（如个别样本系统性低质量）。

---

## 7.导出与保存
需要保存的重要文件为生成的`.html`文件，利用 `cp` 命令保存即可：
```bash
cp qc_results/multiqc_report.html project_summary_QC.html
```
建议保留：
- 各样本 FastQC HTML 报告
- MultiQC 汇总报告
- 修剪前后 FastQC 对比图

> 这些报告可作为分析记录与审稿附录，确保分析可追溯。

---

## 8.总结
- **目的**：FastQC QC 是 RNA-seq 数据分析的第一步，用于评估测序质量、识别接头污染和偏倚。
- **流程**：FastQC → MultiQC 汇总 → （必要时）修剪 → 复检。

- **质量达标标准**：
    - 平均 Q 值 ≥ 30；
    - 接头污染 < 1%；
    - GC 分布符合物种预期；
    - reads 长度稳定。
- **关键输出**：HTML 报告 + 可视化质量曲线。

> 高质量的原始数据是可靠差异分析和富集分析的前提。若 QC 结果异常，应先修复数据再继续下游分析。

---

## 9.常见问题解答（FAQ）
**1.我的 FastQC 报告有 “FAIL”，是不是数据不能用了？**

不一定。FastQC 的阈值较严格。若只在“Per base sequence content”或“Adapter content”模块 FAIL，但其他指标正常，修剪后即可继续。

**2.MultiQC 无法识别 FastQC 结果？**

确认 FastQC 输出目录中包含 .zip 文件，且运行时指定了正确路径：multiqc qc_results/

**3.如何判断是否需要修剪？**

若质量曲线在读段末端明显下滑（<Q20）或接头污染上升，建议修剪。否则可直接进入对齐。

**4.双端（paired-end）样本是否必须同时修剪？**
必须。若只修剪一端，会导致 reads 对齐失配或数量不一致。

**5.修剪会不会影响下游分析？**

合理的修剪（仅去除低质量和接头部分）通常提高比对率；但过度修剪会丢失有效信息，应控制最小长度（如 30 bp 以上）。

**6.常见的接头序列？**

Illumina TruSeq 常见接头：
```bash
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
```
可由 FastQC 的 Overrepresented Sequences 模块检测出。

---