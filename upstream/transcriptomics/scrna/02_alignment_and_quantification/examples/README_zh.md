# 单细胞测序比对与定量（Alignment & Quantification）
在完成原始数据质控（QC）后，下一步是将测序读段（reads）比对到参考基因组并对每个细胞的转录本进行定量。与 bulk RNA-seq 不同，单细胞RNA测序（scRNA-seq）不仅要求准确的比对，还要根据 barcode（细胞标识） 和 UMI（分子标识） 对每个细胞内的转录本去重复计数，最终生成表达矩阵（gene × cell）。

## 1.什么是比对与定量？
### 定义与概念
可参考bulk RNA-seq的概念。

### 比对的核心思路
**scRNA-seq 比对与定量的关键在于识别每条read属于哪个细胞、哪个转录本、哪个分子**：

| 阶段                | 任务                          | 说明             |
| ----------------- | --------------------------- | -------------- |
|  **Barcode提取**     | 从 R1 中读取 Cell barcode 和 UMI | 用于区分不同细胞和转录本分子 |
|  **比对（Alignment）** | 将 R2 读段比对到参考基因组或转录组         | 评估测序准确性与转录区域覆盖 |
|  **UMI去重复与定量**     | 根据 barcode+UMI 合并重复计数       | 生成基因-细胞表达矩阵    |
|  **细胞过滤**          | 移除空油滴、低质量细胞                 | 留下真实细胞用于分析     |


---

## 2.常用工具与依赖
单细胞RNA测序的比对与定量通常通过专用工具完成，这些工具能够识别每条reads的细胞条形码（Cell barcode）和分子标识符（UMI），并进行去重复和表达量统计。
目前主流方案可分为两大类：
| 类别       | 工具                            | 特点                                      | 适用场景                        |                                     |
| -------- | ----------------------------- | --------------------------------------- | --------------------------- | ----------------------------------- |
| **官方方案** | **Cell Ranger（10x Genomics）** | 功能完备、一键式流程、可视化报告                        | 10x Chromium 数据（V2/V3/V3.1） |                                     |
| **开源方案** | **STARsolo / kallisto bustools / salmon alevin / alevin-fry**        |   可定制、轻量、支持多平台                | Drop-seq、inDrops、SMART-seq3、10x自建流程 |

### 核心工具简介
| 工具                             | 主要功能                              | 依赖                      | 输出格式              | 优点          | 适合人群        |          |
| ------------------------------ | --------------------------------- | ----------------------- | ----------------- | ----------- | ----------- | -------- |
| **Cell Ranger**                | 全流程自动化分析（barcode识别、比对、UMI计数、细胞识别） | STAR、Python、C++         | 10x格式矩阵（MTX+TSV）  | 稳定、标准化、报告丰富 | 适合初学者和10x数据 |          |
| **STARsolo**                   | STAR内置模块，支持多种平台                   | STAR、GTF、FASTA          | MTX               | 高效、兼容性强     | 熟悉命令行的科研人员  |          |
| **kallisto bustools**                    |        pseudoalignment + UMI计数                 | kallisto、bustools | MTX     | 极快、内存低      | 大规模数据或教学 |
| **salmon alevin / alevin-fry** | 比对自由定量、UMI校正                      | salmon、rust-based       | MTX / H5AD        | 准确、灵活       | 想自定义流程的用户   |          |
| **seqkit / fastp / samtools**  | 辅助工具，用于FASTQ检查与BAM统计              | 无                       | 文本/表格             | 通用命令行工具     | QC与二次分析     |          |

### 安装
推荐使用 conda 创建虚拟环境后进行安装：
```bash
conda create -n scrna aligners fastqc multiqc # 创建虚拟环境
conda activate scrna # 激活虚拟环境
conda install -c bioconda cellranger star kallisto bustools alevin-fry salmon fastp samtools # 安装依赖
```

类似地，也可以直接下载可执行文件（以Cell Ranger 为例）：
```bash
wget https://cf.10xgenomics.com/releases/cell-exp/cellranger-8.0.0.tar.gz
tar -zxvf cellranger-8.0.0.tar.gz
export PATH=$PATH:/path/to/cellranger-8.0.0
```

### 参考资源
| 工具            | 官方教程/资源                                                                                                                                      |                                                                  |
| ------------- | -------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------- |
| Cell Ranger   | [https://support.10xgenomics.com/single-cell-gene-expression/software](https://support.10xgenomics.com/single-cell-gene-expression/software) |                                                                  |
| STARsolo      | [https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md)             |                                                                  |
| kallisto bustools                                                                                                                                     | [https://www.kallistobus.tools/](https://www.kallistobus.tools/) |
| salmon alevin | [https://salmon.readthedocs.io/en/latest/alevin.html](https://salmon.readthedocs.io/en/latest/alevin.html)                                   |                                                                  |
| alevin-fry    | [https://alevin-fry.readthedocs.io/en/latest/](https://alevin-fry.readthedocs.io/en/latest/)                                                 |                                                                  |

### 工具选择建议
| 需求场景                   | 推荐工具                         | 理由                   |
| ---------------------- | ---------------------------- | -------------------- |
| **标准10x Chromium数据**   | Cell Ranger                  | 自动化、一键报告             |
| **轻量教学或低配计算节点**        | STARsolo / kallisto-bustools | 快速、内存需求低             |
| **大样本批量分析**            | alevin-fry                   | 多线程高效、内存优化           |
| **自定义建库方案（Drop-seq等）** | STARsolo                     | 可灵活定义 barcode/UMI 位置 |
| **可重复性教学展示**           | kallisto-bustools            | 命令简洁，运行快             |


> 在实际科研中，Cell Ranger 是最常用、最稳妥的选择；而在教学、开发和轻量演示中，STARsolo 或 kallisto-bustools 更灵活快速。建议初学者从 Cell Ranger 入手，再逐步掌握 STARsolo 的可配置性。

---

## 3.数据准备与输入
输入文件为下机测序之后的原始文件：
- `R1.fastq.gz`：含有 Cell barcode + UMI
- `R2.fastq.gz`：转录本序列
- `I1.fastq.gz（可选）`：Index序列（文库标签）

目录结构推荐：
```bash
project/
│── raw_data/
│   ├── sample1_R1.fastq.gz
│   ├── sample1_R2.fastq.gz
│── refdata/
│   ├── genome.fa
│   ├── genes.gtf
│── cellranger_output/
│── starsolo_output/
│── qc_results/
│── logs/
```

---

## 4.分析流程
### 4.1 使用 Cell Ranger 进行比对与定量
Cell Ranger 是 10x Genomics 官方提供的分析工具，集成了 barcode识别、比对、UMI去重复、基因计数等功能。输入为原始 FASTQ，输出为标准化的表达矩阵，可直接导入 Seurat 或 Scanpy。

#### 4.1.1 构建参考基因组索引
```bash
cellranger mkref \
  --genome=refdata-gex-GRCh38-2020A \
  --fasta=genome.fa \
  --genes=genes.gtf
```

> 官方提供预构建参考（如人类 GRCh38、小鼠 mm10），也可自行构建。

#### 4.1.2 运行 count 模块（核心步骤）
```bash
cellranger count \
  --id=sample1 \
  --transcriptome=refdata-gex-GRCh38-2020A \
  --fastqs=raw_data/ \
  --sample=sample1 \
  --expect-cells=5000 \
  --localcores=8 \
  --localmem=64
```
输出结果目录及文件：
```bash
cellranger_output/sample1/
│── outs/
│   ├── filtered_feature_bc_matrix/  # 高质量细胞矩阵
│   ├── raw_feature_bc_matrix/       # 未过滤矩阵
│   ├── web_summary.html             # 可视化报告
│   ├── metrics_summary.csv          # 统计指标
│   ├── possorted_genome_bam.bam     # 比对结果
```

> Cell ranger会自动完成：Barcode识别与错误纠正比对（使用 STAR）；UMI计数与过滤；细胞识别与报告生成。

#### 4.1.3 查看报告
使用浏览器打看文件：`web_summary.html`s，其中的报告包含：
- Barcode数与有效细胞数
- 平均reads/细胞
- Mapping比例
- UMI 去重复率
- 空滴油检测（EmptyDrops）结果

### 4.2 使用 STAR solo 进行比对与定量（开源方案）

STARsolo 是 STAR 比对器自带的单细胞RNA-seq模块，轻量高效，支持多平台格式（10x、Drop-seq、Smart-seq3等）。

#### 4.2.1 构建索引
```bash
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir ref_STAR/ \
     --genomeFastaFiles genome.fa \
     --sjdbGTFfile genes.gtf \
     --sjdbOverhang 100
```

#### 4.2.2 执行比对与定量
以10x V3数据为例：
```bash
STAR --runThreadN 8 \
     --genomeDir ref_STAR/ \
     --readFilesIn raw_data/sample1_R2.fastq.gz raw_data/sample1_R1.fastq.gz \
     --readFilesCommand zcat \
     --soloType CB_UMI_Simple \
     --soloCBwhitelist 10x_whitelist.txt \
     --soloCBstart 1 --soloCBlen 16 \
     --soloUMIstart 17 --soloUMIlen 12 \
     --soloFeatures Gene \
     --soloOutFileNames starsolo_output/ sample1
```
输出结果（`starsolo_output/Solo.out/`）：
```bash
│── matrix.mtx         # 稀疏矩阵格式
│── barcodes.tsv       # 细胞ID
│── features.tsv       # 基因信息
│── Summary.csv        # 统计信息
```

> STARsolo 输出的文件格式与 Cell Ranger 完全兼容，可直接加载至 Seurat / Scanpy。

## 5.结果解读
结果分为许多文件：
| 文件                    | 含义          | 下游用途   |
| --------------------- | ----------- | ------ |
| `matrix.mtx`          | 稀疏矩阵（基因×细胞） | 构建表达矩阵 |
| `barcodes.tsv`        | 所有识别的细胞条形码  | 细胞索引   |
| `features.tsv`        | 基因注释信息      | 行标签    |
| `web_summary.html`    | 可视化报告       | 质量评估   |
| `metrics_summary.csv` | 各类统计指标      | 分析记录   |

### 比对质量评估与过滤标准
| 指标                    | 理想范围      | 说明               |
| --------------------- | --------- | ---------------- |
| 比对率（Mapping Rate）     | ≥ 80%     | 表示大部分reads成功比对   |
| UMI去重复率               | ≥ 70%     | 说明分子识别准确         |
| 细胞识别数（Detected cells） | 接近实验预期    | 明显过多说明空油滴        |
| 平均 reads/细胞           | 5万–10万    | 过低影响饱和度          |
| 平均基因数/细胞              | > 1000（人） | 过低提示文库或barcode问题 |
| Mito基因比例              | < 10–15%  | 过高代表细胞应激或死亡      |

## 6.总结

| 阶段    | 工具                        | 主要功能               | 输出     |
| ----- | ------------------------- | ------------------ | ------ |
| QC    | FastQC / MultiQC          | 测序质量检查             | HTML报告 |
| 比对与定量 | CellRanger / STARsolo     | 比对、UMI去重复、细胞识别     | 表达矩阵   |
| 质量评估  | web_summary / Summary.csv | Mapping、UMI统计、细胞过滤 | 指标报告   |

> 比对与定量完成后，我们即可获得标准的**表达矩阵（gene × cell）**，为下游的标准化、降维、聚类与差异分析奠定基础。

---

## 7.常见问题解答（FAQ）
**1.比对率很低（<60%）？**

检查基因组版本、GTF注释与建库试剂是否匹配。若是混合物种实验，需使用联合参考（human + mouse）。

**2.识别到的细胞远少于预期？**

检查Barcode质量（R1 Q分布）、whitelist是否匹配版本。

**3.UMI去重复率过低？**

表示PCR偏倚或Barcode识别不准，可检查接头结构。

**4.空油滴过多？**

Cell Ranger自动用EmptyDrops过滤，可进一步人工调整阈值。

**5.Cell Ranger太重，能否用轻量方案？**

STAR solo或Kallisto/bustools均可，资源需求较低。

---

