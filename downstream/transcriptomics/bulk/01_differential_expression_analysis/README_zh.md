# 差异基因分析
当完成上游分析之后，得到的主要是可用于下游分析的整理好、标准化的矩阵或文件（原始counts/FPKM/TPM文件，**行：基因或转录本；列：样本；内容：原始计数或归一化计数（例如TPM/FPKM）**）。每个样本的均是在一定的处理条件下测定的，所以有一个问题就是处理条件对基因表达有哪些影响，哪些基因表达上调/下调？所以需要进行差异基因分析。


## 1.什么是差异基因分析
差异表达分析（DEA）用于比较不同实验条件或组别下基因的表达水平，识别在特定条件下显著上调或下调的基因。这些基因可能与生物学过程、疾病状态或处理响应相关，是转录组分析中最核心的下游步骤之一。

**应用场景：**
- 不同处理或药物条件下基因表达变化
- 健康与疾病样本比较
- 不同时间点或组织类型的表达动态研究
- 其他对比条件...

> 差异表达分析的核心是统计比较各组基因表达量，判断哪些基因的变化超出了随机波动的范围。常用工具包括R的**DESeq2**、**edgeR**、**limma**等。此外，Python也提供了能够进行差异基因分析的包。

--- 


## 2.工具
### 核心分析工具
- **R / DESeq2**：最常用差异表达分析工具，适合 count 数据，稳定且功能丰富。
- **R / edgeR**：适合小样本或复杂实验设计，基于负二项分布建模。
- **R / limma**：将RNA-seq counts转化为logCPM，适合多条件分析。
- **Python / DESeq2,edgeR,limma via rpy2**：在Python环境调用DESeq2，edgeR,limma。

### 辅助数据处理工具
- **tximport (R)**：将转录本计数汇总为基因计数。
- **biomaRt (R)**：基因注释和ID转换。
- **AnnotationDbi (R)**：本地基因注释数据库操作。
- **pheatmap / ComplexHeatmap (R)**：差异基因可视化。
- **EnhancedVolcano (R)**：绘制火山图。
- **ggplot2 (R)**：通用绘图。
- **clusterProfiler**: 集成包，包含基因注释、ID转换、绘图等。

> 以上的这些工具包都可以通过去观看查看相关的使用手册，或者在R中使用 ``` vigntte(package="xxx")```函数查看、列出和打开R中指定的某个包附带的说明文档。通常比帮助文件（```?function```）更详细。它常包含：使用实例；方法论介绍；代码示例；分析流程展示等。

---

## 3.数据输入
- 基因表达矩阵（gene expression matrix）：经过上游处理得到的counts矩阵（**基因 x 样本**）。
- 样本信息表（colData）：包含样本条件、批次等元信息，例如每个组的处理信息（加药组vs不加药组，处理组vs对照组，或者更多的组别）。

---

## 4.详细流程
### 4.1 DEseq2
#### 步骤1：构建 DESeqDataSet
- **目的**：将count矩阵和样本信息组织成DESeq2可识别的对象。
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
---

## 5.总结
差异基因分析是转录组学研究中最核心的步骤之一，其目标是识别在不同实验条件下显著上调或下调的基因，以揭示潜在的生物学机制或调控网络。整个流程通常包括以下关键环节：
- **1.数据准备**：输入为标准化的基因表达矩阵（如 counts、FPKM、TPM）及样本信息表（包含分组、批次、处理条件等）。
- **2.选择核心的分析工具包**：
    - DESeq2：最常用的差异分析工具，基于负二项分布模型，自动完成归一化、离散度估计与显著性检验，适合大多数 RNA-seq count 数据。
    - edgeR：适合小样本、复杂设计或需要精细建模的情况，同样基于负二项分布。
    - limma：适用于 logCPM 或 microarray 数据，计算效率高，适合多因素分析。
- **3.常规分析步骤**：构建分析对象（如 DESeqDataSet 或 DGEList） -> 数据归一化与低表达过滤 -> 模型拟合与统计检验（估计 Fold Change 与显著性） -> 提取显著差异基因（通常设定阈值 padj < 0.05 且 |log2FC| > 1，也可以根据实验需求自行进行调整） -> 结果展示与解释（火山图：展示差异基因的显著性与表达变化幅度；热图：显示差异基因在样本间的表达模式；MA-plot：展示平均表达量与Fold Change关系；结合生物学知识与功能注释工具（如 clusterProfiler、biomaRt）可进一步进行富集分析与通路解析）。

> **总之，差异基因分析通过系统的统计比较，帮助研究者从大规模转录组数据中筛选关键基因，为后续的生物学验证、功能研究和通路富集分析提供重要依据。**

---

## 6. 常见问题解答（FAQ）
**1.输入数据应该是什么格式？**

差异基因分析通常使用一个原始基因计数矩阵（counts matrix），行为基因、列为样本，内容为每个样本对应的原始计数（未归一化）。
此外还需要一个样本信息表（colData），包含每个样本的分组或条件信息（例如 control / treated）。

> 注意：TPM、FPKM等归一化矩阵不建议直接用于DESeq2或edgeR的差异分析，因为这些方法会自行进行统计归一化。

**2.DESeq2 和 edgeR 有什么区别？我该用哪个？**
| 特点    | DESeq2          | edgeR            |
| ----- | --------------- | ---------------- |
| 统计模型  | 负二项分布           | 负二项分布            |
| 适用数据  | count 数据        | count 数据         |
| 归一化方法 | size factor     | TMM              |
| 优势    | 稳定、适合中大样本量、操作简单 | 灵活、支持复杂设计、小样本表现好 |
| 推荐场景  | 一般差异分析          | 小样本或复杂设计         |

>一般建议：样本量较多、设计简单时用DESeq2；样本量少或需高级模型时用 edgeR。

**3.为什么结果文件里有很多NA或padj为NA？**

这通常由以下原因造成：该基因在大部分样本中表达量太低或全为零；组间方差过大，统计检验无法收敛；样本数过少，统计功效不足。

>解决方法：在分析前过滤掉低表达基因（例如平均 count < 10 的基因）。

**4.padj、pvalue、log2FoldChange 分别是什么意思？**

分别为以下的意思：
- pvalue：单个基因在统计检验中的显著性（未校正）。
- padj（adjusted p-value）：多重检验校正后的p值（通常使用 Benjamini-Hochberg方法），控制假阳性率(FDR)。
- log2FoldChange：基因在两个条件间的表达变化倍数（以2为底的对数，数学表达式为：$\log_2\frac{\text{Gene in condition A}}{\text{Gene in condition B}}$），正值表示上调，负值表示下调。

**5.如何选择差异基因的筛选阈值？**

常用筛选标准为：```padj < 0.05 & abs(log2FoldChange) > 1```，即多重校正后的 p 值小于 0.05 且倍数变化超过 2 倍。当然可以根据研究目的适当放宽或收紧阈值（如 padj<0.01或log2FC>1.5）。

**6.火山图和MA图的区别？**
| 图形类型                         | 横轴                      | 纵轴              | 作用                       |
| ---------------------------- | ----------------------- | --------------- | ------------------------ |
| **火山图 (Volcano Plot)**       | log2FoldChange          | -log10(p-value) | 显示显著性与变化幅度，直观筛选显著差异基因    |
| **MA 图 (Mean-Average Plot)** | 平均表达量 (mean expression) | log2FoldChange  | 显示表达量与差异倍数关系，查看归一化与离散度情况 |

**7.DESeq2 能直接输入FPKM/TPM吗？**

不推荐。DESeq2 和 edgeR 基于原始 count 数据建立模型，FPKM/TPM 已经过归一化，无法正确估计离散度与方差。如果只有FPKM/TPM，可使用limma包进行分析（如voom转换）。

**8.我得到了差异基因列表，接下来该做什么?**

下一步通常包括：
- 功能富集分析：GO、KEGG、Reactome 等（使用 clusterProfiler）。
- 通路可视化：pathview、ggplot2 等。
- 共表达或调控网络分析：WGCNA、Cytoscape。
- 实验验证分析：qPCR、蛋白水平验证等。

> 后续将会一一介绍相应的分析，可以跳转到相应的章节进行学习！
---
