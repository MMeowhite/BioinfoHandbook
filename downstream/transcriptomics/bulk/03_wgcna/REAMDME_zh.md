# WGCNA共表达网络分析

在完成差异基因分析与富集分析后，我们已经能够识别出特定条件下显著变化的基因，并理解它们的功能。但基因之间并非独立工作，它们常常呈现**协同表达（co-expression）**的关系。WGCNA（Weighted Gene Co-expression Network Analysis，加权基因共表达网络分析） 通过表达模式相似性，将基因聚类为模块（modules），并进一步分析这些模块与样本表型的相关性，从而识别潜在的关键模块与核心基因（hub genes）。

## 1.什么是WGCNA？
WGCNA 是一种系统生物学方法，用于从大规模表达数据中识别**基因共表达模块**，并研究模块与表型（如处理组、疾病状态、生理指标）的关系。

核心思想：
- 通过基因间表达模式的相关性构建加权网络。
- 相似表达模式的基因聚为同一个模块（module）。
- 计算每个模块与表型的相关性，筛选与特定性状高度相关的模块。
- 在这些模块中识别出关键调控基因（hub genes）。

---

## 2.常用工具
### 分析与可视化工具（R）

- **WGCNA**：核心加权基因共表达网络分析包，用于构建共表达网络、识别模块、计算模块-表型相关性、提取 hub genes。
- **dynamicTreeCut**：用于动态剪切基因树，识别模块（WGCNA依赖）。
- **igraph**：网络分析与可视化，可对模块内基因网络进行拓扑分析。
- **Cytoscape / RCy3**：将 TOM 网络导出到 Cytoscape，进行交互式可视化与网络分析。
- **pheatmap / ComplexHeatmap**：绘制模块基因热图、样本聚类热图，直观展示模块表达模式。
- **ggplot2**：通用绘图，用于模块特征向量、模块-表型关系和基因-模块关系可视化。
- **clusterProfiler**：可对模块 hub genes 做 GO/KEGG/Reactome 富集分析，支持多种可视化（dotplot、cnetplot、emapplot）。
- **enrichplot**：clusterProfiler 的可视化扩展，绘制模块富集结果图。
- **biomaRt / AnnotationDbi / org.Hs.eg.db**：基因注释与 ID 转换，方便 hub genes 与数据库映射。

--- 

## 3.数据输入与准备
> WGCNA分析一般使用标准化后的表达矩阵（如TPM或vst-transformed counts），而不是原始count。

输入数据：
- **表达矩阵**：行为基因（gene），列为样本（sample）。
- **样本信息表（traits）**：包含每个样本的表型或实验条件（如control/treated、临床信息、连续变量等）。

输入数据要求：
- 建议基因数量：5000–15000（过多计算负担大，可先过滤低方差基因）。
- 样本数量：≥15 个（样本太少会导致网络不稳定）。

---

## 4.分析流程
### 4.1 安装与加载包

```R
# 安装与加载
if(!require(WGCNA)) install.packages("WGCNA")
library(WGCNA)

# 允许多线程
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
```

### 4.2 数据预处理与质量控制
```R
# 读入标准化后的表达矩阵
datExpr <- read.csv("normalized_expression.csv", row.names=1)
datExpr <- as.data.frame(t(datExpr))  # 转置，使行为样本，列为基因

# 检查缺失值与异常样本
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if(!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

# 层次聚类检测异常样本
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
```
> 若有明显异常样本，可根据聚类树手动去除。

### 4.3 软阈值筛选（Soft-threshold selection）
- **目的**：WGCNA 需要选择一个合适的软阈值（power）来构建加权网络，使网络符合无尺度特性（scale-free topology）。
```R
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit", type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, col="red")
abline(h=0.9, col="red")  # 一般选取拟合指数达到0.9的最小power
```
> 选择使拟合指数（$R^2$）约为0.9的最小power值，通常在6–12之间。

### 4.4 构建网络与识别模块
```R
softPower <- 8  # 根据上一步选择的结果设定
adjacency <- adjacency(datExpr, power = softPower)

# 计算Topological Overlap Matrix (TOM)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# 层次聚类识别模块
geneTree <- hclust(as.dist(dissTOM), method = "average")
plot(geneTree, main = "Gene clustering on TOM-based dissimilarity")

# 动态剪切识别模块
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
table(dynamicMods)
```

### 4.5 模块合并与可视化
```R
# 将模块颜色编码
dynamicColors <- labels2colors(dynamicMods)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# 计算模块特征向量并合并相似模块
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes")

merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.25, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs
```

### 4.6 模块与表型的相关性分析
```R
# 读取样本信息
traitData <- read.csv("traits.csv", row.names = 1)
traitData <- as.data.frame(traitData)

# 计算模块与表型的相关性
moduleTraitCor <- cor(mergedMEs, traitData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

# 可视化
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traitData),
               yLabels = names(mergedMEs),
               ySymbols = names(mergedMEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               main = "Module-trait relationships")
```

> 该热图展示了每个模块与样本特征的相关性和显著性。深红（正相关）或深蓝（负相关）模块通常为与表型最强相关的关键模块。

### 4.7 提取关键模块与核心基因（hub genes）
- **目的**：
```R
# 选择与目标性状显著相关的模块
module_of_interest <- "blue"
moduleGenes <- (mergedColors == module_of_interest)

# 计算模块内基因与模块特征值的相关性（Module Membership, MM）
MM <- cor(datExpr[, moduleGenes], mergedMEs[, paste0("ME", module_of_interest)], use="p")

# 计算基因与目标表型的相关性（Gene Significance, GS）
trait <- traitData$Treatment  # 替换为你的表型变量
GS <- cor(datExpr[, moduleGenes], trait, use="p")

# 筛选 hub genes
hub_candidates <- names(sort(abs(MM * GS), decreasing = TRUE))[1:20]
write.csv(hub_candidates, "WGCNA_hub_genes.csv")
```

---

## 5. 结果解读、可视化与导出
### 结果解读
| 分析结果                 | 含义                          |
| -------------------- | --------------------------- |
| **模块（Module）**       | 表达模式相似的基因群体，颜色编码区分          |
| **模块特征向量（ME）**       | 代表该模块整体表达趋势的第一主成分           |
| **Module-Trait 相关性** | 模块与样本表型的显著相关性               |
| **Hub Gene**         | 模块中与ME和表型最强相关的关键基因，潜在核心调控基因 |

### 可视化与导出
- 模块与表型关系热图（`labeledHeatmap()`）
- 模块树状图（`plotDendroAndColors()`）
- 网络可视化（可使用 `TOMplot()` 或导出到 Cytoscape）
```R
# 导出到 Cytoscape
TOM <- TOMsimilarityFromExpr(datExpr, power = softPower)
cyt <- exportNetworkToCytoscape(TOM,
                                edgeFile = "CytoscapeInput-edges.txt",
                                nodeFile = "CytoscapeInput-nodes.txt",
                                weighted = TRUE,
                                threshold = 0.02,
                                nodeNames = colnames(datExpr),
                                altNodeNames = colnames(datExpr),
                                nodeAttr = mergedColors)
```

---

## 6.总结

WGCNA 是差异分析的有力延伸，提供更高层次的系统生物学理解。

分析逻辑： 构建基因共表达网络，识别模块；模块代表了共同表达、可能功能相似的基因群；模块与表型相关性揭示潜在的功能模块；模块内的Hub基因可能是关键调控分子。

优势：不依赖人为设定阈值；能从全局表达模式发现潜在功能模块；有助于挖掘核心基因、标志物与通路。

> 建议进一步与目标性状显著相关模块的Hub基因进行后续功能富集与通路分析（可结合`clusterProfiler`继续深入）。

---

## 7.常见问题解答（FAQ）
**1.WGCNA 可以用原始 count 数据吗？**
- 不推荐直接使用原始 count 数据。WGCNA 基于相关性计算，通常需要 连续型、标准化后的表达矩阵（如 VST 转换的 DESeq2 数据、logCPM 转换的 edgeR 数据或 TPM）。
- 使用原始 count 会导致高度偏态分布，影响网络构建与模块识别。

**2.样本量太少可以做 WGCNA 吗？**
- 样本数过少会导致模块不稳定。一般建议 至少 15 个样本以上。
- 样本太少时，可以考虑先做差异基因筛选，再对高变基因进行 WGCNA，或使用聚合模块方法。

**3. 为什么要选择软阈值（soft-threshold）？**
- WGCNA构建加权网络时需要保证网络近似无尺度拓扑（scale-free topology）。而`pickSoftThreshold()`用于选择合适的power值，使网络符合无尺度特性。
- 选择原则：拟合指数$R^2$ ≥ 0.8–0.9，且平均连通度适中.

**4. 模块颜色和模块特征向量（ME）是什么意思？**
- 模块颜色：表示同一模块基因群（协同表达）。
- 模块特征向量（ME, module eigengene）：模块所有基因表达的第一主成分，代表模块整体表达模式，用于与表型进行相关性分析。

**5.如何选择hub gene？**
- Module Membership（MM）：基因与模块特征向量的相关性。
- Gene Significance（GS）：基因与表型的相关性。
- Hub gene：通常选择MM和GS都较高的基因（如 top 5–20 个），这些基因可能是模块的核心调控基因。

**6.为什么模块与表型相关性很低？**
可能原因：
- 样本量不足，统计功效低；
- 表型与基因表达信号弱或非线性相关；
- 预处理数据不适合构建网络（如过多低表达基因或未归一化数据）。

**7. 是否需要去除异常样本？**
- 需要，异常样本会影响网络构建和模块识别。
- 解决方法：可通过样本聚类（`hclust()`）检查异常样本并决定是否去除。

**8. WGCNA 可以分析多条件或连续型表型吗？**
- 可以。表型可以是 二分类、连续型或多类别。
- 对于多因素实验设计，可在表型矩阵中添加多个列，模块-表型相关性将逐列计算。

**9. 模块可以进一步进行功能富集分析吗？**

完全可以。对hub genes或模块基因进行GO/KEGG/Reactome富集分析，可以发现模块潜在的生物学功能。推荐使用`clusterProfiler`、`ReactomePA`等包进行可视化和注释。

**10. WGCNA 分析需要多大计算资源？**
基因数量较多（>10,000）或样本量较大时，计算TOM和聚类比较耗时。可以使用`allowWGCNAThreads()`多线程加速。同时对低方差基因先进行过滤，减少计算量。

---