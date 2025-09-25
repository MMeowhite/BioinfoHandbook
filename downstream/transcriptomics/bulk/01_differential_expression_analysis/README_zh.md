# 差异基因分析
当完成上游分析之后，得到的主要是可用于下游分析的整理好、标准化的矩阵或文件（原始counts/FPKM/TPM文件，行：基因或转录本；列：样本；内容：原始计数或归一化计数（例如TPM/FPKM））。每个样本的均是在一定的处理条件下测定的，所以有一个问题就是处理条件对基因表达有哪些影响，哪些基因表达上调/下调？所以需要进行差异基因分析。


## 1.什么是差异基因分析
差异表达分析（DEA）用于比较不同实验条件或组别下基因的表达水平，识别在特定条件下显著上调或下调的基因。这些基因可能与生物学过程、疾病状态或处理响应相关，是转录组分析中最核心的下游步骤之一。

**应用场景：**
- 不同处理或药物条件下基因表达变化
- 健康与疾病样本比较
- 不同时间点或组织类型的表达动态研究

差异表达分析的核心是统计比较各组基因表达量，判断哪些基因的变化超出了随机波动的范围。常用工具包括R的**DESeq2**、**edgeR**、**limma**等。


## 2.工具
### 核心分析工具：
- **R / DESeq2**：最常用差异表达分析工具，适合 count 数据，稳定且功能丰富。
- **R / edgeR**：适合小样本或复杂实验设计，基于负二项分布建模。
- **R / limma**：将RNA-seq counts转化为logCPM，适合多条件分析。
- **Python / DESeq2,edgeR,limma via rpy2**：在Python环境调用DESeq2，edgeR,limma。

### 辅助数据处理工具：
- **tximport (R)**：将转录本计数汇总为基因计数。
- **biomaRt (R)**：基因注释和ID转换。
- **AnnotationDbi (R)**：本地基因注释数据库操作。
- **pheatmap / ComplexHeatmap (R)**：差异基因可视化。
- **EnhancedVolcano (R)**：绘制火山图。
- **ggplot2 (R)**：通用绘图。
- **clusterProfiler**: 集成包，包含基因注释、ID转换、绘图等。


## 3.数据输入
- 基因表达矩阵（gene expression matrix）：经过上游处理得到的counts矩阵（基因 x 样本）。
- 样本信息表（colData）：包含样本条件、批次等元信息，例如每个组的处理信息（加药组vs不加药组，处理组vs对照组，或者更多的组别）。


## 4.详细流程
### 4.1 DEseq2
#### 步骤1：构建 DESeqDataSet
- **目的**：将 count 矩阵和样本信息组织成DESeq2可识别的对象。
- **示例代码**：
```R
# counts：原始基因计数矩阵，不需要提前归一化
# colData：样本信息表，至少包含 condition 列表示实验组
# design：设计公式，表示比较的条件。多因素实验可以写成 ~ batch + condition
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ condition)
```
- **注意**：避免使用全为零的基因，DESeq2会自动过滤低表达基因，但也可以手动过滤（例如count<10）

#### 步骤2：差异表达分析
- **目的**：计算基因在不同条件下的统计显著性和Fold Change
- **示例代码**：
```R
dds <- DESeq(dds)
res <- results(dds)
```
- **说明**：DESeq() 会自动进行归一化（size factor）、离散度估计和 Wald 检验；results() 可以指定比较组，例如res <- results(dds, contrast=c("condition","treated","control"))

#### 步骤3：可视化结果
- MA-plot：显示基因表达Fold Change与平均表达量的关系
- 火山图 (Volcano plot)：显示 Fold Change 与统计显著性（-log10 p-value）
- 热图 (Heatmap)：展示差异基因在样本间的表达模式

#### 步骤4：导出差异基因列表
- **目的**：将显著差异基因保存为表格，便于下游分析
- **示例代码**：
```R
diff_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1) #取多重检验校正后的p-value<0.05以及log2FoldChange的绝对值>1的基因
write.csv(as.data.frame(diff_genes), file="diff_genes.csv")
```
- **说明**：可以根据需要调整阈值或排序输出；也可以只取上/下调的基因


### 4.2 edgeR
#### 步骤1：准备数据
- 输入：counts matrix (genes × samples)，样本信息表 `group` 
- 示例代码：
```R
library(edgeR)

dge <- DGEList(counts = counts, group = colData$condition)
dge <- calcNormFactors(dge)  # TMM 归一化
```

#### 步骤2：数据过滤
- **目的**：移除低表达基因，提高统计功效
- **示例代码**：
```R
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]
```

#### 步骤3：构建设计矩阵
- **目的**：构建可以被edgeR设别的对象
- **示例代码**：
```R
design <- model.matrix(~0 + group, data=colData)
colnames(design) <- levels(colData$condition)
dge <- estimateDisp(dge, design)
```

#### 步骤4：差异分析
- **目的**：计算基因在不同条件下的统计显著性和Fold Change
- **示例代码**：
```R
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, contrast=c(-1,1))  # 对比条件：treated vs control
topTags(qlf)
```

#### 步骤5：可视化结果
- 还是MA-plot，火山图，热图等
- **示例代码**：
```R
plotMD(qlf)                     # MA-plot
hist(qlf$table$PValue)          # P-value 分布
```

#### 步骤6：导出差异基因
- **示例代码**：
```R
diff_genes <- topTags(qlf, n=Inf)
write.csv(diff_genes, file="edgeR_diff_genes.csv")
```
