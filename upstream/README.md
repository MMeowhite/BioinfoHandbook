# 上游分析
上游分析（upstream analysis）一般指从原始测序数据（通常是 FASTQ）到结构化、可定量的生物学矩阵或表格的全过程。它主要关注**数据质量控制、对齐、定量、特征提取，而不是生物学解释**。

> 下游分析（downstream analysis）则是从这些结果出发，进行**差异分析、聚类、注释、通路富集、可视化等**。

## 1.什么是上游分析？
上游分析的核心思想：**上游分析 = 把“原始测序数据（序列）”变成“可分析的数字表格/文件”的整个技术过程，即把测序仪输出的原始文件（一般是 FASTQ）处理成下游可用的标准化结果（例如基因×样本的 counts 矩阵、VCF 变异文件、peak 文件、质谱的蛋白丰度表等）。**

为什么上游分析这么重要？ 
- 质量问题（如测序坏、接头、污染）如果不在上游修掉，会在下游产生假阳性/假阴性；
- 不同上游策略（比对工具、注释版本）会严重影响可重复性和结果解读；
- 好的上游输出能让下游分析更稳定、解释更可信。

---

## 2.常见类型的生信数据与对应的上游分析

| 数据类型                                | 典型原始数据        | 上游分析主要步骤                         | 核心输出                    | 常用工具                                        |
| ----------------------------------- | ------------- | -------------------------------- | ----------------------- | ------------------------------------------- |
| **DNA测序 (WGS/WES)**                 | FASTQ         | QC → 比对 → 去重 → 变异检测 → 注释         | VCF（变异文件）               | BWA, GATK, samtools, bcftools, FastQC       |
| **转录组 (bulk RNA-seq)**              | FASTQ         | QC → 修剪 → 比对 → 定量 → 归一化          | count矩阵 / TPM表          | FastQC, STAR, HISAT2, featureCounts, Salmon |
| **单细胞转录组 (scRNA-seq)**              | FASTQ (R1/R2) | QC → barcode解析 → 比对 → 去重 → 定量    | count矩阵 (cell × gene)   | Cell Ranger, STARsolo, Alevin               |
| **表观组 (ATAC-seq)**                  | FASTQ         | QC → 比对 → 去重 → peak calling      | peak文件 (BED/narrowPeak) | Bowtie2, MACS2, deepTools                   |
| **单细胞ATAC (scATAC-seq)**            | FASTQ         | QC → barcode解析 → 比对 → peak矩阵     | cell × peak矩阵           | Cell Ranger ATAC, ArchR                     |
| **甲基化测序 (WGBS / RRBS)**             | FASTQ         | QC → 比对（双链识别）→ 甲基化位点调用           | CpG methylation表        | Bismark, bsmap, MethyKit                    |
| **小RNA测序 (miRNA-seq)**              | FASTQ         | QC → 去接头 → 比对到miRBase → 计数       | miRNA表达矩阵               | miRDeep2, bowtie, sRNAbench                 |
| **ChIP-seq**                        | FASTQ         | QC → 比对 → 去重 → peak calling → 富集 | peak表 + 信号轨迹            | Bowtie2, MACS2, HOMER                       |
| **空间转录组 (Spatial transcriptomics)** | FASTQ + 图像    | QC → 比对 → 空间定位 → 定量              | 空间矩阵 + 组织图              | Space Ranger, STutility                     |
| **宏基因组 / 16S**                      | FASTQ         | QC → 组装/分类 → OTU/ASV定量           | OTU/ASV表 + 物种注释         | QIIME2, Kraken2, MetaPhlAn                  |
| **蛋白质组 (LC-MS/MS)**                 | mzML, RAW     | 峰检测 → 定量 → 归一化                   | 蛋白丰度表                   | MaxQuant, Proteome Discoverer               |
| **代谢组**                             | mzML          | 峰提取 → 对齐 → 定量 → 注释               | 代谢物丰度表                  | XCMS, MZmine, MS-DIAL                       |

## 3.上游分析包括哪些“步骤”？（通用版）
上游分析拆成 6 个泛化的步骤，几乎适用于任何组学类型（RNA、DNA、ATAC、蛋白质组等）：

### 3.1 原始数据接收与组织（Data Acquisition & Organization）
- **目的**：确保原始文件完整、命名规范、路径清晰，为后续批量处理做好准备。

#### 常见文件类型
| 数据类型                     | 原始文件格式           | 示例文件名                                    |
| ------------------------ | ---------------- | ---------------------------------------- |
| Bulk RNA-seq / scRNA-seq | FASTQ / FASTQ.GZ | sample1_R1.fastq.gz, sample1_R2.fastq.gz |
| ATAC-seq / ChIP-seq      | FASTQ / FASTQ.GZ | atac_R1.fastq.gz                         |
| WGS / WES                | FASTQ / FASTQ.GZ | WGS_001_R1.fastq.gz                      |
| 蛋白质组（质谱）                 | RAW / MZML       | sample1.raw                              |
| 甲基化测序                    | FASTQ            | sample_R1.fastq.gz                       |

常见操作：
- 从测序中心下载 .fastq.gz 或 .tar.gz 文件；
- 使用 md5sum 校验文件完整性；
- 重命名样本文件（统一命名规则，如 sample_condition_replicate_R1.fastq.gz）；

组织目录结构：
```bash
project/
├── raw_data/          # 原始FASTQ文件
├── qc_results/        # 质控报告
├── trimmed_data/      # 修剪后FASTQ
├── alignment/         # 比对结果BAM
├── quantification/    # 定量结果
└── logs/              # 日志与统计文件
```

> **小技巧**：一个项目必须有一个**metadata表（sample_sheet.csv）**，记录样本名、分组信息、测序批次、read数、性别、组织来源等，方便后续批量处理与追溯。


### 3.2 原始质量控制（QC）
- **目的**：评估测序数据的整体质量、检测污染、判断是否需要修剪或重新测序。
- **典型工具**：
    - FastQC：单样本质控（生成 HTML 报告）
    - MultiQC：整合所有样本 QC 报告
    - fastp：结合 QC 与修剪功能（可生成交互式报告）

#### 常用检查指标
| 模块                  | 说明      | 异常信号       | 处理建议  |
| ------------------- | ------- | ---------- | ----- |
| Per base quality    | 每碱基质量分布 | 末端下降（Q<20） | 修剪末端  |
| GC content          | GC含量分布  | 双峰/偏移      | 检查污染  |
| Adapter content     | 接头污染    | 明显上升       | 去接头   |
| Overrepresented seq | 重复序列    | 接头/rRNA污染  | 过滤或修剪 |

输出文件：
- `*_fastqc.html / .zip`
- `multiqc_report.html`（汇总报告）

> 结论判断：若大部分碱基在 Q30 以上、接头污染 <1%，即可进入下一步。

### 3.3 预处理（Trimming / Filtering）
- **目的**：去除低质量碱基、接头序列、污染片段，保证后续比对准确。
- **常用工具**：
    - TrimGalore（Cutadapt 封装）
    - fastp（速度快、功能全）
    - Trimmomatic（经典老牌）

典型命令：
```bash
fastp -i sample_R1.fastq.gz -I sample_R2.fastq.gz \
      -o trimmed_R1.fastq.gz -O trimmed_R2.fastq.gz \
      --detect_adapter_for_pe --thread 8 --html fastp_report.html
```

注意事项：
- 单细胞测序数据：不要修剪 R1（包含 cell barcode + UMI），只可处理 R2。
- 设置最小长度（`--length_required 30`），避免过度修剪导致 reads 丢失。

输出文件：
- 修剪后 FASTQ（通常更小）
- 修剪报告（`fastp.html / .json`）

### 3.4 对比或伪对比
- **目的**：将 reads 定位到参考基因组或转录组，确定每条 reads 的来源（哪个基因/区域）。

| 分析类型                | 目的               | 常用工具                                | 输出格式           |
| ------------------- | ---------------- | ----------------------------------- | -------------- |
| Bulk RNA-seq        | 转录组比对或伪比对        | STAR / HISAT2 / Salmon / kallisto   | BAM / quant.sf |
| scRNA-seq           | barcode-aware 比对 | Cell Ranger / STARsolo / alevin-fry | matrix.mtx     |
| ATAC-seq / ChIP-seq | 基因组比对            | Bowtie2 / BWA / STAR                | BAM            |
| DNA-seq（WGS/WES）    | 基因组比对            | BWA / minimap2                      | BAM            |
| 甲基化测序               | CpG定位比对          | Bismark                             | BAM            |

示例命令（RNA-seq, STAR）：
```bash
STAR --runThreadN 8 \
     --genomeDir ref_index \
     --readFilesIn sample_R1.trim.gz sample_R2.trim.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix alignment/sample_ \
     --outSAMtype BAM SortedByCoordinate
```
伪比对（Salmon / kallisto）：
```bash
salmon quant -i ref_index -l A \
      -1 sample_R1.trim.gz -2 sample_R2.trim.gz \
      -p 8 -o quant/sample
```
> **比对率参考值：RNA-seq ≥70%；DNA-seq ≥90%；ATAC-seq 约 50–70%**。


### 3.5 去重复 / UMI处理 / 计数（Dedup / UMI handling / Quantification）
- **目的**：消除 PCR 扩增偏倚，统计每个基因或区域的有效 reads 数量。

不同数据类型的处理方法：
| 数据类型         | 操作内容                         | 常用工具                                |
| ------------ | ---------------------------- | ----------------------------------- |
| Bulk RNA-seq | 从 BAM 定量到基因                  | featureCounts / htseq-count         |
| scRNA-seq    | UMI 聚合、去重复、生成 cell × gene 矩阵 | Cell Ranger / alevin-fry / STARsolo |
| ATAC-seq     | 去重复后 peak calling            | Picard MarkDuplicates / MACS2       |
| DNA-seq      | 标记重复、变异识别                    | Picard / GATK                       |
| 蛋白质组         | 光谱定量与蛋白归一                    | MaxQuant / FragPipe                 |

示例（RNA 定量）：
```bash
featureCounts -T 8 -p -B -C -a gencode.gtf \
              -o counts.txt alignment/*.bam
```

输出：
- 基因 × 样本 counts 表
- 单细胞矩阵（matrix.mtx / barcodes.tsv / features.tsv）
- 去重复 BAM / metrics 文件

> 小贴士：在单细胞数据中，UMI 是去重的关键，确保 R1 中的 UMI 正确解析。

### 3.6 格式化、汇总与归档（Output）
- **目的**：将分析结果标准化保存，方便下游使用与结果复现。

输出类型（按数据类型）：
| 数据类型         | 上游输出                                     | 下游输入（示例）          |
| ------------ | ---------------------------------------- | ----------------- |
| Bulk RNA-seq | gene_counts.txt                          | DESeq2 / edgeR    |
| scRNA-seq    | matrix.mtx + barcodes.tsv + features.tsv | Seurat / Scanpy   |
| ATAC-seq     | BAM + peaks.bed                          | DiffBind / ArchR  |
| DNA-seq      | BAM + VCF                                | Mutect2 / ANNOVAR |
| 蛋白质组         | proteinGroups.txt                        | Perseus / MSstats |


推荐组织结构：
```bash
project/
├── raw_data/ 
├── qc_results/
├── trimmed_data/
├── alignment/
├── quantification/
├── reports/       # MultiQC、日志汇总
└── summary/       # 最终结果（counts、matrix、VCF等）
```
归档建议：
- 保留所有日志文件（便于复现）
- 记录软件版本（建议使用 conda env export > env.yml）
- 对大型文件（BAM/FASTQ）可压缩存储或上传至云端（如 SRA）

### 3.6 小结
| 步骤         | 目的         | 主要输出           | 关键工具                        |
| ---------- | ---------- | -------------- | --------------------------- |
| 原始数据接收 | 样本组织与核对    | FASTQ          | —                           |
| QC     | 检查测序质量     | HTML报告         | FastQC / MultiQC            |
| 预处理    | 去接头与低质量序列  | 修剪FASTQ        | fastp / TrimGalore          |
| 对比     | 映射reads到参考 | BAM / quant.sf | STAR / BWA / Salmon         |
| 计数     | 基因或UMI定量   | counts矩阵       | featureCounts / Cell Ranger |
| 汇总     | 输出标准化结果    | 表格/矩阵/VCF      | MultiQC / 自定义脚本             |

> 总结：上游分析是从“原始序列”到“数字化表达矩阵”的桥梁，是所有生物信息学分析的起点。

## 不同组学上游分析的共性
虽然研究对象不同，但上游分析通常包含这些共同环节：

| 阶段                     | 功能                           | 说明                            |
| ---------------------- | ---------------------------- | ----------------------------- |
| **1. 原始数据 QC**         | 检查测序质量、接头污染                  | FastQC / MultiQC              |
| **2. 修剪与过滤**           | 去除低质量 reads、接头               | fastp / Trimmomatic           |
| **3. 比对或伪比对**          | 将 reads 映射到参考基因组 / 转录组 / 数据库 | BWA, STAR, kallisto, Bowtie2  |
| **4. 去重复 / 去嵌合 / 去污染** | 去除PCR重复或宿主污染                 | samtools rmdup / DecontaMiner |
| **5. 特征提取 / 定量**       | 计算每个基因、peak、位点的信号强度          | featureCounts, MACS2          |
| **6. 格式化与汇总**          | 统一输出为可下游分析的矩阵                | count矩阵, VCF, BED, 表达表等       |

## 上游分析常见生态系统（Pipeline框架）
| 框架                             | 适用领域                    | 说明                                |
| ------------------------------ | ----------------------- | --------------------------------- |
| **nf-core**                    | 多组学自动化分析                | 提供 WGS / RNA / scRNA / ATAC 等标准流程 |
| **Cell Ranger / Space Ranger** | 10x系列（单细胞、空间）           | 10x官方上游分析套件                       |
| **bcbio-nextgen**              | RNA / DNA / ChIP / ATAC | 全自动化 pipeline 平台                  |
| **snakemake / nextflow**       | 通用框架                    | 可自定义任意上游分析工作流                     |

## 上游分析需要准备/掌握的技能与资源
上游分析较为复杂，但是一般利用现成的工具跑Pipeline即可。
- 基础命令行（Linux）：文件操作、管道、脚本编写
- 常用生信工具：FastQC、samtools、bwa/STAR、salmon、featureCounts、GATK、MACS2、Cell Ranger（看数据类型）
- 工作流框架（可选但推荐）：Nextflow、Snakemake（便于自动化与复现）
- 环境管理：Conda、Docker/Singularity（保证工具版本稳定）
- 计算资源：WGS 需要大量 CPU/内存/磁盘；单细胞和 bulk RNA 也需要中等资源。建议先用小数据做测试。

## 总结

> 上游分析是“从序列到数字”的过程；下游分析是“从数字到生物学意义”的过程。

---

