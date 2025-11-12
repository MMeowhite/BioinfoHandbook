# SCENIC转录因子调控网络分析
SCENIC（Single-Cell rEgulatory Network Inference and Clustering）是一种用于单细胞 RNA-seq 数据的**转录因子（TF）调控网络推断与活性评分分析**方法。 它可以从基因表达矩阵中识别出 活跃的转录调控模块（regulons），并计算每个细胞的TF活性，从而揭示细胞类型特异性或状态依赖的调控机制。

## 1.什么是SCENIC转录因子调控网络分析？
SCENIC是一种针对单细胞 RNA-seq 数据的转录因子调控网络分析方法。它的核心目标是：**从单细胞基因表达矩阵中，推断出哪些转录因子（TF）在不同细胞群体中处于活跃状态，并控制哪些下游靶基因（target genes）**。通俗理解：普通的 scRNA-seq 只能告诉你“哪些基因表达了”，而 SCENIC 能告诉你“为什么这些基因表达了”【是由哪些转录因子（如 MYC、GATA3、NF-κB 等）调控的】。
### 常见应用场景：
- **细胞类型特异性转录因子识别**：判断哪些 TF 主导不同 cluster 的表达程序。
- **发育/分化过程中的调控动态**：可结合 RNA velocity 轨迹展示 TF 活性变化趋势。
- **肿瘤细胞与免疫微环境相互作用**：分析免疫调控 TF（如 STAT、NF-κB）在不同细胞中的差异。
- **与 bulk 数据或 ChIP-seq 结果整合**：验证调控关系的实验支持。

### 分析思路与原理：
SCENIC 主要分为三个核心步骤：
| 步骤             | 工具                     | 主要功能                                 |
| -------------- | ---------------------- | ------------------------------------ |
| **网络推断**     | **GENIE3 / GRNBoost2** | 基于基因共表达构建 TF → target 预测网络           |
| **motif 分析** | **RcisTarget**         | 从预测的 TF–target 关系中筛选具有保守 motif 的目标基因 |
| **活性评分**     | **AUCell**             | 计算每个细胞中各 regulon 的活性得分（AUC 值）        |

最终输出的：
- Regulon 列表（TF 及其靶基因集合）
- AUC 矩阵（每个细胞的 regulon 活性）

--- 

## 2.常用工具与安装
| 工具包                    | 语言       | 功能                |
| ---------------------- | -------- | ----------------- |
| **SCENIC**             | R        | 官方 R 实现           |
| **pySCENIC**           | Python   | 更快，适合大规模数据        |
| **GENIE3 / GRNBoost2** | R/Python | 基因调控网络推断          |
| **RcisTarget**         | R        | motif 富集分析        |
| **AUCell**             | R        | regulon 活性计算      |
| **SCENICplus**         | Python   | 支持多组学（ATAC+RNA）扩展 |

---

## 3.数据输入与准备
输入数据：
- 单细胞表达矩阵（如 Seurat 对象或 `.loom` 文件）
- 基因名需为 gene symbol
- 需要对应物种的 motif 数据库（`hg38`, `mm10` 等）

示例文件准备：
```R
# 从 Seurat 对象导出表达矩阵
library(Seurat)
expr_mat <- GetAssayData(seurat_obj, slot = "data")
write.csv(expr_mat, file = "expr_matrix.csv")
```

## 4.分析流程
### 4.1 R 分析流程
#### 4.1.1 加载库与数据导入
```R
library(SCENIC)
library(GENIE3)
library(RcisTarget)
library(AUCell)

exprMat <- read.csv("expr_matrix.csv", row.names = 1)
cellInfo <- read.csv("cell_metadata.csv", row.names = 1)
```

#### 4.1.2 构建 co-expression 网络
```R
genesKept <- geneFiltering(exprMat, minCountsPerGene = 3*0.01*ncol(exprMat))
exprMat_filtered <- exprMat[genesKept, ]
runGenie3(exprMat_filtered, nTrees = 1000)
```

#### 4.1.3 motif 富集与 regulon 构建
```R
scenicOptions <- initializeScenic(org = "hgnc", dbDir = "cisTarget_databases")
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)
```

#### 4.1.4 输出 regulon 活性矩阵
```R
aucell <- loadInt(scenicOptions, "aucell_regulonAUC")
saveRDS(aucell, file = "SCENIC_regulon_AUC.rds")
```

### 4.2 Python 分析流程
#### 4.2.1 基因共表达网络构建
```bash
pyscenic grn expr_matrix.csv TF_list.txt -o adj.tsv --num_workers 8
```

#### 4.2.2 motif 富集分析
```bash
pyscenic ctx adj.tsv cisTarget_db/*.feather --annotations_fname motifs.tbl \
  --expression_mtx expr_matrix.csv --mode "dask_multiprocessing" \
  -o regulons.csv --num_workers 8
```

#### 4.2.3 计算 regulon 活性
```bash
pyscenic aucell expr_matrix.csv regulons.csv -o auc_mtx.csv --num_workers 8
```

## 5.结果解读与可视化
SCENIC 的核心结果包括 转录因子调控模块（regulons）、每个细胞的调控活性（AUC 值） 以及 调控网络结构。这些结果帮助我们理解不同细胞群体背后的转录调控逻辑。

### 5.1 Regulon 活性（AUC）分析与可视化
AUC 矩阵（Regulon activity matrix） 反映每个 regulon（转录因子调控模块）在每个细胞中的活性。
常见可视化方式包括热图、UMAP 映射、箱线图等。
实例代码：
```R
library(SCENIC)
library(ComplexHeatmap)
library(ggplot2)

# Regulon AUC 矩阵
regulonAUC <- getAUC(scenicOptions)
regulonAUC_mat <- t(as.data.frame(regulonAUC))

# 热图展示高变 regulon 活性
top_regulons <- names(sort(apply(regulonAUC_mat, 2, sd), decreasing = TRUE))[1:20]
Heatmap(regulonAUC_mat[, top_regulons],
        name = "Regulon Activity (AUC)",
        show_row_names = FALSE,
        cluster_columns = TRUE)
```
解读：
- 不同细胞群体在某些 regulon 上呈现特异性活性；
- 活性较高的 TF 通常为该群体的关键调控因子；
- 可进一步结合聚类结果或伪时间分析，观察 TF 活性随状态变化的动态趋势。

### 5.2 Regulon 活性在 UMAP 上的空间分布
可将 regulon 活性值映射到单细胞 UMAP，以展示特定 TF 在不同细胞群体中的活性差异。
```R
# 选择感兴趣的 TF
tf_interest <- "GATA3"

# 将 regulon 活性加入 Seurat 对象
seurat_obj[["GATA3_AUC"]] <- regulonAUC_mat[tf_interest, colnames(seurat_obj)]

# 可视化 TF 活性分布
FeaturePlot(seurat_obj, features = "GATA3_AUC", cols = c("lightgrey", "red"))
```
解读：
- 高活性区域通常对应 TF 主导的细胞亚群；
- **若 TF 在轨迹分析中随伪时间上升，可能为关键命运决定因子**。

### 5.3 Regulon 与细胞类型关系
利用平均 AUC 值展示不同细胞类型或聚类的 regulon 活性差异。
```R
# 计算每个 cluster 的平均 AUC 值
cluster_auc <- data.frame(
  cluster = seurat_obj$seurat_clusters,
  t(regulonAUC_mat)
) %>%
  group_by(cluster) %>%
  summarise_all(mean)

# 可视化（条形图或热图）
pheatmap(as.matrix(cluster_auc[, -1]),
         cluster_rows = TRUE, cluster_cols = TRUE,
         annotation_row = data.frame(Cluster = cluster_auc$cluster))
```
解读：
- 某些 TF 在特定 cluster 中显著上调，提示它们是该细胞群体的核心调控因子；
- 可与 marker gene 结果进行交叉验证。

### 5.4 转录调控网络可视化
SCENIC 输出的 regulon 网络可导入 Cytoscape 进行可视化，展示 TF 与其靶基因的调控关系。
```R
# 导出网络到 Cytoscape
export2cytoscape(scenicOptions, minGenes = 5, weightThreshold = 0.1)
```
> 在 Cytoscape 中可进行：网络拓扑分析（度中心性、介数中心性），筛选 hub TF（高连接度的关键转录因子），模块化展示（Regulon-based subnetworks）

解读：
- 高度连接的 TF 是潜在的全局调控因子；
- 同模块内 TF–target 关系可能构成细胞类型特异的调控程序。

### 5.5 其他整合分析
SCENIC可以与其他分析整合一起揭示更加深层的生物调控机制。
| 整合对象                          | 整合目的                                 |
| ----------------------------- | ------------------------------------ |
| **轨迹分析（Monocle / Slingshot）** | 观察 TF 活性在发育伪时间上的变化，识别关键命运调控因子        |
| **CellChat / CellPhoneDB**    | 将 TF 活性与信号通路富集结果关联，揭示转录调控与通讯网络的层级关系  |
| **差异基因与富集分析**                 | 对高活性 TF 的靶基因集进行 GO/KEGG 富集，探究调控的功能方向 |

---

## 6.总结
SCENIC 是单细胞 RNA-seq 分析中用于**转录因子调控网络推断与活性评估**的关键工具。它通过整合基因共表达、motif 调控信息和细胞层面的活性评分，能够揭示细胞类型或状态背后的转录调控机制。

分析逻辑：
- 共表达网络构建：利用 GENIE3 / GRNBoost2 推断 TF 与下游靶基因的关系；
- motif 验证：通过 RcisTarget 识别真实可能结合的基因集合（Regulon）；
- 细胞层面活性计算：用 AUCell 评估每个 regulon 在单个细胞中的活性（AUC 值）；
- 可视化与解释：基于 regulon 活性矩阵（AUC matrix）重新聚类、绘制热图或与细胞类型关联。

主要输出：
| 结果类型           | 含义                     |
| -------------- | ---------------------- |
| **Regulon**    | 由转录因子及其靶基因组成的调控单元      |
| **AUC Matrix** | 每个细胞的 regulon 活性得分矩阵   |
| **TF 网络图**     | 展示 TF–target 调控关系的网络结构 |
| **细胞调控特征图**    | 各细胞中 TF 活性分布，揭示群体差异    |

SCENIC 的优势在于：能够在单细胞水平揭示调控层面的异质性；输出具有生物学解释力的 TF–target 网络；结果可与聚类结果、轨迹分析、CellChat 通讯结果等进行整合；对发现关键转录因子、细胞命运决定因子特别有价值。

> 总之，SCENIC 不仅能看见基因表达，更能**解释谁在驱动这些表达**。
它是连接单细胞转录组分析与调控机制研究的重要桥梁。


## 7.常见问题解答（FAQ）
**1.SCENIC 与 WGCNA 有什么区别？**

WGCNA 是基于 bulk RNA-seq 的共表达模块，SCENIC 是基于单细胞数据的 TF 调控网络。WGCNA 强调“基因共表达”，而 SCENIC 强调“调控方向性”。

**2.R 版还是 Python 版更好？**
主要取决于你的数据集大小。小数据集（<5k 细胞）：R 版即可。大数据集：推荐 `pySCENIC`，速度快 10–20 倍。

**3.SCENIC 分析需要哪些数据库？**

主要需要以下的数据库：
- motif 数据库 (hg38, mm10 等 feather 文件)
- TF 注释表（可从官方 GitHub 获取）
- cisTarget motif annotation 文件

**4.AUCell 得分代表什么？**

表示某个 TF 对应的 regulon 在单个细胞中的活跃程度（AUC 值越高 → 活性越强）。

**5.结果与 cluster 不对应？**

需要注意的是 SCENIC 关注调控活性，不一定与表达相同；可用 UMAP 或 PCA 对 regulon AUC 矩阵重新聚类。

---