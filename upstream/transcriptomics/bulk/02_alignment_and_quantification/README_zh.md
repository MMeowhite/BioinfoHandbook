# 比对与定量（Alignment & Quantification）
在完成原始数据质控（QC）与必要的修剪后，RNA-seq 分析的下一步是将测序读段（reads）比对到参考基因组或转录组，并统计每个基因或转录本的表达量。
该步骤直接决定了下游差异分析的精度与可靠性，是从“原始序列”到“可量化基因表达矩阵”的关键环节。

## 1.什么是比对与定量？
### 比对
比对是指**将测序得到的短序列（reads）映射到参考基因组或转录组，确定每个 read 的来源位置**。
RNA-seq 特有的“跨内含子拼接”特性要求使用专门的spliced aligner（剪接感知比对工具），如：
- HISAT2
- STAR
- TopHat2

### 定量
定量是指**根据比对结果（或伪比对结果），统计每个基因或转录本的表达量**。
可以采用两类方法：
- **基于比对的定量**：例如 featureCounts、HTSeq-count；
- **伪比对（quasi-mapping）**：如 Salmon、Kallisto（更快，准确率高，教学友好）。

> 比对负责“定位”，定量负责“计数”；两者是 RNA-seq 数据转化为表达矩阵的关键步骤。

## 2.常用工具
| 任务     | 工具                              | 特点                 |
| ------ | ------------------------------- | ------------------ |
| 比对到基因组 | **STAR / HISAT2**               | 支持剪接感知；输出 BAM 文件   |
| 基因计数   | **featureCounts / HTSeq-count** | 基于注释的基因计数          |
| 伪比对+定量 | **Salmon / Kallisto**           | 无需生成 BAM；速度快，资源消耗低 |
| 结果整合   | **tximport (R)**                | 汇总转录本计数为基因层级表达矩阵   |

---

## 3.数据输入与准备
输入文件：
- 修剪后的 FASTQ 文件（`.fasq.gz`）
- 参考基因组或转录组序列（FASTA）
- 基因注释文件（GTF / GFF3）

输出文件：
- 对齐结果（BAM / SAM）
- 表达定量文件（如 `quant.sf`）
- 汇总的基因表达矩阵（用于差异分析）

同样的，需要进行文件树管理：
```bash
project/
│── trimmed_data/         # 修剪后的FASTQ
│── reference/            # 基因组或转录组文件
│── alignment_results/    # 比对输出(BAM)
│── quant_results/        # 定量结果(Salmon等)
│── counts_matrix/        # 表达矩阵（汇总后）
```

---

## 4.分析流程
### 4.1 参考索引构建
基因组比对，HISAT2为例：
```bash
# 构建索引
hisat2-build reference/genome.fa reference/genome_index
```
伪对比，以Salmon为例：
```bash
salmon index -t reference/transcripts.fa -i reference/salmon_index
```

> 转录本序列（`transcripts.fa`）通常可从Ensembl或Gencode下载

### 4.2 比对reads（Alignments）
#### 4.2.1 HISAT
```bash
hisat2 -x reference/genome_index \
       -1 trimmed_data/sample1_R1.fastq.gz \
       -2 trimmed_data/sample1_R2.fastq.gz \
       -S alignment_results/sample1.sam \
       --rna-strandness RF -p 8

# SAM 转 BAM + 排序
samtools view -bS alignment_results/sample1.sam | samtools sort -o alignment_results/sample1.sorted.bam
samtools index alignment_results/sample1.sorted.bam
```

#### 4.2.2 STAR
```bash
# 先生成索引（一次性操作）
STAR --runThreadN 8 --runMode genomeGenerate \
     --genomeDir reference/star_index \
     --genomeFastaFiles reference/genome.fa \
     --sjdbGTFfile reference/annotation.gtf \
     --sjdbOverhang 100

# 运行比对
STAR --runThreadN 8 --genomeDir reference/star_index \
     --readFilesIn trimmed_data/sample1_R1.fastq.gz trimmed_data/sample1_R2.fastq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix alignment_results/sample1_ \
     --outSAMtype BAM SortedByCoordinate
```

> 输出的 `.bam` 文件可用于 featureCounts 计数或可视化（IGV）

### 4.3 基因定量（Quatification）
#### 4.3.1 基于 BAM 的 featureCounts
```bash
featureCounts -T 8 -p -t exon -g gene_id \
  -a reference/annotation.gtf \
  -o counts_matrix/gene_counts.txt \
  alignment_results/*.sorted.bam
```
输出为基因x样本的计数矩阵，可直接进入 `DESeq2`

#### 4.3.2 伪比对
```bash
salmon quant -i reference/salmon_index \
  -l A \
  -1 trimmed_data/sample1_R1.fastq.gz \
  -2 trimmed_data/sample1_R2.fastq.gz \
  -p 8 \
  -o quant_results/sample1
```
输出文件：
```bash
quant_results/sample1/quant.sf
```
该文件包含：
| 列名       | 含义               |
| -------- | ---------------- |
| Name     | 转录本 ID           |
| Length   | 转录本长度            |
| TPM      | 每百万转录本数          |
| NumReads | 分配到该转录本的 reads 数 |

### 4.4 汇总为基因表达矩阵（R / tximport）
```bash
library(tximport)
library(readr)
samples <- read.table("samples.txt", header=TRUE)   # 包含样本名与路径
files <- file.path("quant_results", samples$SampleID, "quant.sf")
names(files) <- samples$SampleID

# 转录本到基因的映射表（tx2gene）
tx2gene <- read.csv("reference/tx2gene.csv")  # 两列：transcript_id, gene_id

txi <- tximport(files, type="salmon", tx2gene=tx2gene)
countData <- txi$counts
write.csv(countData, file="counts_matrix/gene_counts_matrix.csv")
```

> `tximport` 自动汇总转录本层面数据为基因层面表达矩阵，供 `DESeq2` 或 `edgeR` 使用。

--- 

## 5.结果解读与质量评估
### 5.1 比对质量评估
| 指标                           | 理想范围 | 含义                  |
| ---------------------------- | ---- | ------------------- |
| **总比对率（Overall alignment rate）** | >80% | reads 成功比对到参考基因组的比例 |
| **唯一比对率（Uniquely mapped）**       | >70% | 唯一定位的 reads 占比      |
| **多重比对率（Multiple mapped）**       | <10% | 多位置匹配的 reads 占比     |
| **未比对（Unmapped）**                | <10% | 未能比对的 reads         |

生成统计（HISAT2 示例输出）：
```bash
...
Overall alignment rate = 91.6%
```

MultiQC汇总对齐报告：
```bash
multiqc alignment_results/ -o alignment_results/
```
可快速比较各样本的比对率、reads数与质量。

### 5.2 定量结果检查
- TPM 分布：样本间总体分布应相近；
- 高表达基因是否符合预期（如 housekeeping genes）；
- 低表达基因比例是否正常（极低计数过多可能为测序深度不足）。

### 5.3 下游衔接
完成比对与定量并且质量评估合格之后，可以进入：
- 差异表达分析（DESeq2 / edgeR）
- 可视化（PCA、火山图）
- 功能富集分析（GO/KEGG/GSEA）
- ...

> 这是 RNA-seq 分析从“reads”到“表达值”的转折点。

---

## 6.总结
- 目的：将 RNA-seq 原始 reads 转换为可量化的基因表达矩阵。
- 方法选择：
    - 比对法：STAR / HISAT2 + featureCounts；
    - 伪比对法：Salmon / Kallisto（推荐教学用）。
- 关键步骤：构建参考索引 -> 比对 reads -> 定量基因表达 -> 汇总为基因矩阵 -> 检查比对率与分布。
- 输出结果：
    - 对齐文件（.bam）；
    - 定量文件（quant.sf / counts.txt）；
    - 表达矩阵（counts_matrix.csv）。

> **可靠的比对与定量是后续差异分析与功能富集的基础**。高质量的表达矩阵决定了统计分析的灵敏度与可信度。

---

## 7.常见问题解答（FAQ）
**1.STAR / HISAT2 的比对率太低怎么办？**
- 检查参考基因组与注释是否匹配（Ensembl vs Gencode）；
- 检查样本是否污染或物种错误；
- 检查修剪是否过度导致 read 太短。

**2.Salmon 输出的 quant.sf 文件能直接用于 DESeq2 吗？**

不行，需要通过 tximport 汇总为基因层级矩阵。

**3.featureCounts 的 -t 和 -g 参数代表什么？**
- `-t`：计数的注释类型（通常为 exon）
- `-g`：按哪种注释分组（通常为 gene_id）

**4.单端测序是否也可用 Salmon？**

可以，只需指定单端参数：
```bash
salmon quant -i index -r sample.fastq.gz -l A -o quant/
```

**5.Salmon / Kallisto 与 HISAT2 + featureCounts 哪个更好？**

伪比对（Salmon/Kallisto）速度快、占用低，教学友好；比对法（HISAT2/STAR）更传统，可生成 BAM 供可视化与其他用途。值得注意的是两者在表达定量精度上两者差异极小，所以一般不影响后续的分析结果。

**6.tx2gene 文件从哪来？**

从注释文件（GTF）中提取 transcript_id 与 gene_id 对照关系；可使用R包 `GenomicFeatures`：
```bash
txdb <- makeTxDbFromGFF("reference/annotation.gtf")
k <- keys(txdb, keytype="TXNAME")
tx2gene <- select(txdb, keys=k, keytype="TXNAME", columns="GENEID")
```

---