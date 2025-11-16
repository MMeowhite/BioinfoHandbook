# Snakemake 自动化流程：RNA-seq 上游分析

## 1.项目目录结构
```bash
project/
│── raw_data/            # 原始 FASTQ
│── trimmed_data/        # 修剪后的 FASTQ
│── qc_results/          # FastQC / MultiQC 输出
│── reference/           # 基因组/转录组 + 注释文件
│── alignment_results/   # HISAT2/STAR 输出 BAM
│── quant_results/       # Salmon/Kallisto 输出
│── counts_matrix/       # 汇总表达矩阵
│── config.yaml          # 配置文件
│── Snakefile            # Snakemake 主文件
```

---

## 2.snakemake示例
首先写 `.yaml`配置文件：
```yaml
samples:
  - sample1
  - sample2

# 原始 FASTQ 文件
raw_dir: "raw_data"
trimmed_dir: "trimmed_data"

# 输出目录
qc_dir: "qc_results"
bam_dir: "alignment_results"
quant_dir: "quant_results"
counts_dir: "counts_matrix"

# 基因组 / 转录组
genome_fasta: "reference/genome.fa"
transcriptome_fasta: "reference/transcripts.fa"
annotation_gtf: "reference/annotation.gtf"
tx2gene: "reference/tx2gene.csv"

# HISAT2 / Salmon 索引
hisat2_index: "reference/genome_index"
salmon_index: "reference/salmon_index"

# 线程数
threads: 8
```

然后就是写snakemake的 `.py`脚本：

```python
import os
configfile: "config.yaml"

SAMPLES = config["samples"]
THREADS = config["threads"]

RAW_DIR = config["raw_dir"]
TRIM_DIR = config["trimmed_dir"]
QC_DIR = config["qc_dir"]
BAM_DIR = config["bam_dir"]
QUANT_DIR = config["quant_dir"]
COUNTS_DIR = config["counts_dir"]

HISAT2_INDEX = config["hisat2_index"]
SALMON_INDEX = config["salmon_index"]
ANNOTATION = config["annotation_gtf"]
TX2GENE = config["tx2gene"]

# 默认目标：生成基因表达矩阵
rule all:
    input:
        os.path.join(COUNTS_DIR, "gene_counts_matrix.csv")

# ---------------------------
# 1. FastQC 原始数据质控
rule fastqc_raw:
    input:
        r1=lambda wildcards: f"{RAW_DIR}/{wildcards.sample}_R1.fastq.gz",
        r2=lambda wildcards: f"{RAW_DIR}/{wildcards.sample}_R2.fastq.gz"
    output:
        html1=os.path.join(QC_DIR, "{sample}_R1_fastqc.html"),
        html2=os.path.join(QC_DIR, "{sample}_R2_fastqc.html")
    threads: 2
    shell:
        """
        fastqc -t {threads} -o {QC_DIR} {input.r1} {input.r2}
        """

# ---------------------------
# 2. Trimmomatic 修剪
rule trimmomatic:
    input:
        r1=lambda wildcards: f"{RAW_DIR}/{wildcards.sample}_R1.fastq.gz",
        r2=lambda wildcards: f"{RAW_DIR}/{wildcards.sample}_R2.fastq.gz"
    output:
        r1_paired=os.path.join(TRIM_DIR, "{sample}_R1_paired.fastq.gz"),
        r1_unpaired=os.path.join(TRIM_DIR, "{sample}_R1_unpaired.fastq.gz"),
        r2_paired=os.path.join(TRIM_DIR, "{sample}_R2_paired.fastq.gz"),
        r2_unpaired=os.path.join(TRIM_DIR, "{sample}_R2_unpaired.fastq.gz")
    threads: 4
    shell:
        """
        trimmomatic PE -threads {threads} \
        {input.r1} {input.r2} \
        {output.r1_paired} {output.r1_unpaired} \
        {output.r2_paired} {output.r2_unpaired} \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """

# ---------------------------
# 3. FastQC 修剪后质控
rule fastqc_trimmed:
    input:
        r1=lambda wildcards: f"{TRIM_DIR}/{wildcards.sample}_R1_paired.fastq.gz",
        r2=lambda wildcards: f"{TRIM_DIR}/{wildcards.sample}_R2_paired.fastq.gz"
    output:
        html1=os.path.join(QC_DIR, "{sample}_R1_trimmed_fastqc.html"),
        html2=os.path.join(QC_DIR, "{sample}_R2_trimmed_fastqc.html")
    threads: 2
    shell:
        """
        fastqc -t {threads} -o {QC_DIR} {input.r1} {input.r2}
        """

# ---------------------------
# 4. HISAT2 比对 + BAM 排序
rule hisat2_align:
    input:
        r1=lambda wildcards: f"{TRIM_DIR}/{wildcards.sample}_R1_paired.fastq.gz",
        r2=lambda wildcards: f"{TRIM_DIR}/{wildcards.sample}_R2_paired.fastq.gz",
        index=expand(HISAT2_INDEX + ".1.ht2", 0)
    output:
        bam=os.path.join(BAM_DIR, "{sample}.sorted.bam")
    threads: THREADS
    shell:
        """
        hisat2 -x {HISAT2_INDEX} -1 {input.r1} -2 {input.r2} \
        -S {wildcards.sample}.sam -p {threads} --rna-strandness RF
        samtools view -bS {wildcards.sample}.sam | samtools sort -o {output.bam}
        samtools index {output.bam}
        rm {wildcards.sample}.sam
        """

# ---------------------------
# 5. Salmon 伪比对与定量
rule salmon_quant:
    input:
        r1=lambda wildcards: f"{TRIM_DIR}/{wildcards.sample}_R1_paired.fastq.gz",
        r2=lambda wildcards: f"{TRIM_DIR}/{wildcards.sample}_R2_paired.fastq.gz",
        index=SALMON_INDEX
    output:
        directory=os.path.join(QUANT_DIR, "{sample}")
    threads: THREADS
    shell:
        """
        salmon quant -i {input.index} -l A \
        -1 {input.r1} -2 {input.r2} -p {threads} \
        -o {output.directory}
        """

# ---------------------------
# 6. 汇总为基因表达矩阵
rule tximport_matrix:
    input:
        quant_files=expand(os.path.join(QUANT_DIR, "{sample}/quant.sf"), sample=SAMPLES),
        tx2gene=TX2GENE
    output:
        csv=os.path.join(COUNTS_DIR, "gene_counts_matrix.csv")
    run:
        import rpy2.robjects as ro
        ro.r('library(tximport)')
        samples = SAMPLES
        files = [f"{QUANT_DIR}/{s}/quant.sf" for s in samples]
        files_r = "c(" + ",".join([f'"{f}"' for f in files]) + ")"
        names_r = "c(" + ",".join([f'"{s}"' for s in samples]) + ")"
        tx2gene = TX2GENE
        ro.r(f"""
        files <- {files_r}
        names(files) <- {names_r}
        tx2gene <- read.csv("{tx2gene}")
        txi <- tximport(files, type="salmon", tx2gene=tx2gene)
        write.csv(txi$counts, file="{output.csv}")
        """)
```

然后可以运行以下的命令
- 运行全部流程：
```bash
snakemake -j 8
```

- 只运行某一步（如HISAT2比对）：
```bash
snakemake -j 8 alignment_results/sample1.sorted.bam
```

- 检查文件依赖与DAG：
```bash
snakemake --dag | dot -Tpdf > dag.pdf
```

- 输出结果：
    - **qc_results/**：原始与修剪后的 FastQC 报告
    - **trimmed_data/**：修剪后的 FASTQ
    - **alignment_results/*.sorted.bam**：HISAT2 比对结果
    - **quant_results/*/quant.sf**：Salmon 转录本计数
    - **counts_matrix/gene_counts_matrix.csv**：汇总基因表达矩阵

---