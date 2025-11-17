# Snakemake 自动化流程：单细胞 RNA-seq（scRNA-seq）上游分析

## 1. 项目目录结构
```bash
project/
│── raw_data/                # 原始 FASTQ
│── qc_results/              # fastqc / multiqc 输出
│── cellranger_results/      # cellranger count 输出
│── merged_matrix/           # 合并后的表达矩阵 (h5)
│── reference/               # Cell Ranger reference
│── config.yaml              # 配置文件
│── Snakefile                # Snakemake 主文件
```

---

## 2. Snakemake 示例
### 2.1 config.yaml 配置文件
```yaml
# 样本名称
samples:
  - sample1
  - sample2

# 原始数据
raw_dir: "raw_data"
qc_dir: "qc_results"

# Cell Ranger 输出
cr_count_dir: "cellranger_results"

# 合并输出
merged_dir: "merged_matrix"

# Cell Ranger reference
cellranger_ref: "reference/10x_ref"

# 技术设置（10x）
chemistry: "auto"

# threading
threads: 16
```

### 2.2 Snakemake 工作流程（Snakefile）
```python
import os
configfile: "config.yaml"

SAMPLES = config["samples"]
THREADS = config["threads"]

RAW_DIR = config["raw_dir"]
QC_DIR = config["qc_dir"]
CR_DIR = config["cr_count_dir"]
MERGED_DIR = config["merged_dir"]

CR_REF = config["cellranger_ref"]
CHEM = config["chemistry"]

# 默认目标：合并矩阵
rule all:
    input:
        os.path.join(MERGED_DIR, "merged_matrix.h5")


# -------------------------------------------------
# 1. FastQC 原始数据质控
rule fastqc_raw:
    input:
        fq=lambda wc: expand(f"{RAW_DIR}/{wc.sample}_S1_L001_R{r}_001.fastq.gz", r=[1,2])
    output:
        expand(os.path.join(QC_DIR, "{sample}_R{r}_fastqc.html"), r=[1,2])
    threads: 4
    shell:
        """
        fastqc -t {threads} -o {QC_DIR} {input}
        """


# -------------------------------------------------
# 2. Cell Ranger count（比对 + UMI 去重 + 细胞 barcode 呼叫）
rule cellranger_count:
    input:
        fq1=lambda wc: f"{RAW_DIR}/{wc.sample}_S1_L001_R1_001.fastq.gz",
        fq2=lambda wc: f"{RAW_DIR}/{wc.sample}_S1_L001_R2_001.fastq.gz"
    output:
        directory(os.path.join(CR_DIR, "{sample}"))
    threads: THREADS
    shell:
        """
        cellranger count \
            --id={wildcards.sample} \
            --transcriptome={CR_REF} \
            --fastqs={RAW_DIR} \
            --sample={wildcards.sample} \
            --chemistry={CHEM} \
            --localcores={threads} \
            --localmem=64

        mv {wildcards.sample} {CR_DIR}/
        """


# -------------------------------------------------
# 3. 使用 Cell Ranger aggr 合并多个样本
rule cellranger_aggr:
    input:
        aggr_csv="aggr.csv"
    output:
        directory(os.path.join(MERGED_DIR, "aggregation"))
    threads: THREADS
    shell:
        """
        cellranger aggr \
            --id=aggregation \
            --csv={input.aggr_csv} \
            --normalize=mapped \
            --localcores={threads} \
            --localmem=64

        mv aggregation {MERGED_DIR}/
        """


# -------------------------------------------------
# 4. 提取最终 h5 表达矩阵
rule extract_matrix_h5:
    input:
        directory(os.path.join(MERGED_DIR, "aggregation"))
    output:
        h5=os.path.join(MERGED_DIR, "merged_matrix.h5")
    shell:
        """
        cp {input}/outs/filtered_feature_bc_matrix.h5 {output.h5}
        """
```

---

## 3.aggr.csv 文件格式（用于 cellranger aggr）
你需要手动准备一个 aggr.csv：
```csv
sample_id,molecule_h5
sample1,cellranger_results/sample1/outs/molecule_info.h5
sample2,cellranger_results/sample2/outs/molecule_info.h5
```
放在根目录即可。

---

## 4.运行流程
### 4.1 运行全部流程
```bash
snakemake -j 16
```

### 4.2 只跑某一步（例如 cellranger count）
```bash
snakemake cellranger_results/sample1
```

### 4.3 查看依赖图 DAG
```bash
snakemake --dag | dot -Tpdf > dag_scRNA.pdf
```

---

## 5. 输出结果说明
| 目录                                 | 内容                                  |
| ---------------------------------- | ----------------------------------- |
| **qc_results/**                    | FASTQ 质控报告 (FastQC)                 |
| **cellranger_results/**            | 每个 sample 的比对结果 + 表达矩阵              |
| **merged_matrix/merged_matrix.h5** | 最终合并后的单细胞表达矩阵，可直接用于 Seurat / Scanpy |
| **merged_matrix/aggregation/**     | 完整的 Cell Ranger aggr 输出             |
