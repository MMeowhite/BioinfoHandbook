# 单细胞测序原始数据质控（Quality Control, QC）
在单细胞RNA测序（scRNA-seq）分析流程中，原始数据质控（QC）是确保下游比对、定量和细胞识别可靠性的首要步骤。单细胞数据较bulk RNA-seq更复杂，**不仅包含转录本序列信息，还在读段中嵌入了细胞条形码（Cell barcode）与分子标识符（UMI）**。因此，QC 既要评估测序质量，也要确认 barcode 与 UMI 的完整性和均匀性。

## 1.什么是单细胞测序QC？
在单细胞RNA测序（scRNA-seq）中，每个细胞通过唯一的barcode标识，其转录本的每个分子带有UMI用于去重复计数。
质控的主要目标是：
- 检查测序质量是否合格（Q分布、GC含量、接头污染等）；
- 评估barcode和UMI读段的质量是否足以支持准确的细胞识别；
- 检查数据完整性与测序深度是否符合预期。

与bulk RNA-seq不同，scRNA-seq的FASTQ文件通常分为：
- **R1**：存放细胞条形码（Cell barcode）和UMI序列（通常长度较短）；
- **R2**：存放转录本（mRNA）序列；
- （部分平台还会有 I1 / I2：Index序列，用于文库索引）。

---

## 2.常用工具与依赖
QC部分使用的工具与Bulk RNA-seq类似，如下所示：
| 工具                                | 功能             | 说明                  |
| --------------------------------- | -------------- | ------------------- |
| **FastQC**                        | 基础测序质量分析       | 适用于R1、R2、I1等FASTQ文件 |
| **MultiQC**                       | 汇总多个样本FastQC报告 | 支持多通道数据比较           |
| **seqkit / zcat**                 | 快速查看FASTQ头部信息  | 检查barcode与read结构    |
| **fastp / cutadapt / TrimGalore** | （可选）接头与低质量过滤   | 通常不建议过度修剪单细胞数据      |

安装方法也是可以通过 conda 快速安装以及直接下载可执行文件。

---

## 3.数据输入与准备
### 常见测序数据结构
以10x Genomics为例：
| 文件类型                             | 内容说明                   | 典型长度     | 示例                         |
| --------------------------| ---------------------- | -------- | -------------------------- |
| `sample_S1_L001_R1_001.fastq.gz` | **Cell barcode + UMI**     | 26bp     | TCGTGCGTGTTGACTGAGTAAGTCTG |
| `sample_S1_L001_R2_001.fastq.gz` | **转录本读段（Transcript read）** | 90–150bp | ACGTGGCGTATGCTAGCTGG...    |
| `sample_S1_L001_I1_001.fastq.gz` | **索引序列（Index read）**       | 8bp      | ATCGTTAG                   |

> **注意：R1 的碱基质量直接影响后续 CellRanger/STARsolo 对 barcode 的识别率，是 scRNA QC 的关键指标**。

### 输入文件格式
单细胞测序的原始文件为 `.fastq` 或 `.fastq.gz`，命名通常包含 lane 与 read 信息：
```bash
@A00123:56:H23HLDSXX:1:1101:1100:1000 1:N:0:CGTACTAG
AGCTGATCGTACGTCGATCGATCG...
+
IIIIIIIIIIIIIIIIIIIIIIII...
```

### 项目文件树建议
```bash
project/
│── raw_data/
│   ├── sample1_S1_L001_R1_001.fastq.gz
│   ├── sample1_S1_L001_R2_001.fastq.gz
│── qc_results/
│── trimmed_data/       # 若进行修剪
│── logs/
```

> 建议使用标准命名方式（R1为barcode/UMI，R2为转录本），方便自动化分析。

---

## 4.分析流程
### 4.1 单样本质控（FastQC）
单样本执行 FastQC，输出每个文件的质量报告：
```bash
mkdir -p qc_results
fastqc raw_data/sample1_S1_L001_R1_001.fastq.gz -o qc_results/
fastqc raw_data/sample1_S1_L001_R2_001.fastq.gz -o qc_results/
```
输出结果为：
- `sample1_S1_L001_R1_001_fastqc.html`
- `sample1_S1_L001_R2_001_fastqc.html`
- `.zip` 文件（详细数据）

### 4.2 批量质控（多个样本）
```bash
fastqc raw_data/*.fastq.gz -o qc_results/ -t 8
```
参数说明：
- `-o` 输出路径；
- `-t` 并行线程数（建议4-8）。

### 4.3 多样本汇总报告（MultiQC）
整合所有样本的 FastQC 报告，便于统一评估：
```bash
multiqc qc_results/ -o qc_results/
```
输出：
- `multiqc_report.html`
- `multiqc_data/`（解析结果）

---

## 5.结果解读
### 5.1 FastQC 报告重点（R1 vs R2）
| 模块                        | R1（barcode/UMI） | R2（transcript）  | 异常含义             |
| ------------------------- | --------------- | --------------- | ---------------- |
| Per base quality          | 整体质量应均匀         | 前后段略降属正常        | 若严重下降，影响条形码识别或比对 |
| Sequence length           | 固定（26bp左右）      | 固定（90–150bp）    | 不一致说明拼接问题        |
| Adapter content           | 应接近0            | 可略有升高           | 可能存在接头污染         |
| GC content                | 均匀              | 应符合转录组分布        | 异常峰提示污染          |
| Overrepresented sequences | 条形码序列重复正常       | 若高比例重复说明 PCR 偏倚 | 过多重复可影响UMI计数     |

### 5.2 MultiQC 报告
MultiQC 提供：
- 所有样本的平均Q值对比；
- 每read类型（R1/R2）的质量趋势；
- GC分布对比；
- Adapter污染比例汇总。

> 若不同lane或样本质量差异明显，建议重点复查低质量样本，必要时重新测序。

### 5.3 QC 结果评估标准
| 指标        | 理想值                    | 说明             |
| --------- | ---------------------- | -------------- |
| 平均Q值      | ≥ Q30                  | 质量良好           |
| Adapter污染 | < 1%                   | 过高则需修剪         |
| GC含量分布    | 单峰、平滑                  | 双峰提示混合或污染      |
| R1长度      | 固定                     | 异常波动需复查        |
| Reads数    | 与目标细胞数匹配（约5万 reads/细胞） | 测序深度不均会影响后续饱和度 |

> 简单来说，若 R1（barcode）质量良好，R2 质量中等偏上，即可进入 CellRanger 或 STARsolo 比对步骤。

---

## 6.下游处理
### 6.1 是否需要修剪？
单细胞数据通常不建议进行过度修剪，因为：
- R1 读段中含barcode和UMI，不应改变；
- 修剪可能破坏barcode结构导致比对失败；
- CellRanger/STARsolo会自动识别接头，不需额外cutadapt处理。

> 例外：若FastQC报告中R2接头污染严重，可仅对R2执行轻度修剪。

修剪示例代码：
```bash
fastp -i raw_data/sample1_R2.fastq.gz \
      -o trimmed_data/sample1_R2.trimmed.fastq.gz \
      --detect_adapter_for_pe --thread 4 \
      --html trimmed_data/sample1_fastp.html
```

### 6.2 重现 QC 修剪结果
若执行了修剪，建议重新运行 FastQC 以验证结果
```bash
fastqc trimmed_data/*.fastq.gz -o qc_results/trimmed_qc/
multiqc qc_results/trimmed_qc/ -o qc_results/
```

---

## 7.导出与保存
建议保存以下的文档：
- 每个文件的 .html FastQC 报告；
- MultiQC 汇总报告；
- （如有修剪）修剪前后FastQC结果；
- 质控日志（记录命令、时间、版本）。

参考命令：
```bash
cp qc_results/multiqc_report.html reports/QC_summary_scRNA.html
```

---

## 8.总结

| 关键点  | 内容                      |
| ---- | ----------------------- |
| QC目标 | 检查测序质量、条形码完整性、污染情况      |
| 核心工具 | FastQC、MultiQC          |
| 输入文件 | R1（barcode+UMI）、R2（转录本） |
| 关注重点 | R1质量、GC分布、adapter污染     |
| 不建议  | 随意修剪R1；过度过滤             |
| 输出   | FastQC & MultiQC 报告     |

> 高质量的R1和R2读段是正确识别细胞、准确量化转录本的基础。QC合格后，可进入下游 比对与定量（第2章） 步骤。

---

## 9.常见问题解答（FAQ）
**1.R1 很短（26bp），FastQC 报告很多 “FAIL”，怎么办？**

属于正常现象。R1 包含 barcode 与 UMI，不反映转录本质量；FastQC 只按read长度标准评分，会误判。

**2.需要对R1进行修剪吗？**

不需要！R1序列结构固定（barcode+UMI），修剪会破坏条形码识别。

**3.MultiQC 显示 GC 双峰？**

多数情况是混合物种或rRNA污染，需确认文库制备过程。

**4.R1 与 R2 的 read 数不一致？**

检查是否合并错误或文件丢失。双端测序必须成对存在。

**5.单细胞样本 FastQC 结果普遍比 bulk 差？**

属正常。单细胞建库量少、reads短，Q30略低不影响分析。

**6.修剪后 CellRanger 识别的细胞数变少？**

原因是 barcode 被破坏。**建议只修剪 R2 或完全不修剪**。

---