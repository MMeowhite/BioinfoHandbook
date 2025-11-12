# 肿瘤免疫微环境分析（TME）
肿瘤免疫微环境（Tumor Microenvironment, TME） 是指肿瘤组织中除肿瘤细胞本身以外的所有非恶性成分及其相互作用的综合系统。
它包括免疫细胞（T 细胞、B 细胞、NK 细胞、巨噬细胞、树突状细胞等）、成纤维细胞、血管内皮细胞、基质、细胞因子以及代谢环境等。TME 的复杂免疫状态对肿瘤发生、发展及免疫治疗反应具有重要影响。而单细胞 RNA-seq 提供了刻画 TME 细胞组成、相互作用与功能状态的理想工具。

## 1.什么是肿瘤免疫微环境分析？
从功能性的角度来说，根据大量的研究文章显示肿瘤免疫微环境不仅影响肿瘤的生长、侵袭与转移，还决定了肿瘤对免疫治疗（如 PD-1/PD-L1 抑制剂）的响应。
因此，解析 TME 的细胞组成、功能状态与通讯网络，是理解肿瘤免疫逃逸机制和指导免疫治疗的关键。

因此分析以下的内容至关重要：
| 分析层面      | 研究重点                               |
| ----------- | ---------------------------------- |
| **细胞组成**  | 区分免疫细胞类型与比例（T、B、Mφ、NK、DC等）         |
| **细胞状态**  | 判断细胞活化、耗竭、抑制或极化状态                  |
| **信号通路**  | 分析免疫相关通路（IFN、TGF-β、NF-κB等）         |
| **免疫检查点** | 研究 PD-1、CTLA4、LAG3 等基因表达模式         |
| **细胞通讯**  | 揭示肿瘤–免疫、免疫–免疫之间的配体–受体交互            |
| **调控机制**  | 识别调控免疫细胞状态的关键转录因子（如 IRF、STAT、NFAT） |

所以分析目标主要分为以下的几类：
- 区分免疫细胞类型与亚群（T细胞、B细胞、巨噬细胞等）
- 描述免疫细胞功能状态（活化、耗竭、抑制）
- 分析免疫检查点（PD-1、CTLA4、LAG3）表达
- 解析细胞间通讯（重点是免疫细胞与肿瘤细胞）
- 识别关键调控因子（转录因子、信号通路）

---

## 2.常用工具与数据库
常用的分析工具（R）：
| 工具                         | 功能简介                   |
| -------------------------- | ---------------------- |
| **Seurat**                 | 单细胞数据整合、聚类与免疫细胞注释基础    |
| **SingleR / scCATCH**      | 免疫细胞自动注释               |
| **CellChat / CellPhoneDB** | 细胞通讯分析                 |
| **SCENIC / DoRothEA**      | 转录因子调控网络推断             |
| **GSVA / PROGENy**         | 通路与信号活性分析              |
| **CIBERSORTx / EPIC**      | 免疫细胞比例估算（适合 bulk 数据对照） |

常用的免疫数据库：
| 数据库                                 | 内容             |
| ----------------------------------- | -------------- |
| **CellMarker / PanglaoDB**          | 各类免疫细胞标记基因     |
| **ImmPort / MSigDB (immune sets)**  | 免疫相关基因集        |
| **TISIDB / TIMER / TCIA**           | 免疫浸润与免疫治疗相关数据库 |
| **Reactome / KEGG Immune Pathways** | 免疫信号通路注释       |

---

## 3.数据输入与准备
数据来源与适用分析类型：
| 数据类型             | 可进行的 TME 分析               |
| ---------------- | ------------------------- |
| **bulk RNA-seq** | 免疫细胞比例估算（CIBERSORTx、EPIC） |
| **scRNA-seq**    | 免疫亚群识别、功能状态、细胞通讯、调控网络     |
| **空间转录组（ST）**    | 免疫细胞空间分布与肿瘤细胞交互           |
| **多组学整合**        | 结合突变、甲基化、代谢组，揭示调控机制       |

## 4.分析流程
常规的分析思路比较简单：
单细胞表达矩阵 -> 免疫细胞识别与注释 -> 免疫亚群细分与 marker 基因提取 -> 免疫功能与通路活性分析（GSVA / PROGENy）-> 免疫检查点表达分析 ->  细胞通讯分析（CellChat）-> 转录因子调控网络（SCENIC）

> 当然也可以加入个性化定制分析，记住生信分析永远是需要围绕你想要揭示的问题进行解析。

### 4.1 免疫细胞识别与提取
```R
# 提取免疫相关 cluster
immune_clusters <- subset(sc_data, idents = c("T_cells", "B_cells", "Macrophages"))

# 重新聚类
immune_clusters <- FindClusters(immune_clusters, resolution = 0.5)
DimPlot(immune_clusters, reduction = "umap", label = TRUE)
```

### 4.2 亚群注释与 marker 鉴定
```R
# 自动注释
library(SingleR)
pred <- SingleR(test = as.SingleCellExperiment(immune_clusters),
                ref = HumanPrimaryCellAtlasData(),
                labels = ref$label.main)
immune_clusters$SingleR_labels <- pred$labels
```

### 4.3 功能通路活性分析（GSVA）
```R
library(GSVA)
gsva_res <- gsva(as.matrix(GetAssayData(immune_clusters, slot = "data")),
                 immune_pathways,
                 method = "ssgsea")

# 可视化结果
heatmap(gsva_res)
```

### 4.4 免疫检查点基因表达
```R
checkpoint_genes <- c("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2")
FeaturePlot(immune_clusters, features = checkpoint_genes, cols = c("lightgrey", "red"))
```

### 4.5 细胞通讯分析（CellChat）
```R
library(CellChat)
cellchat <- createCellChat(object = sc_data, group.by = "celltype")
cellchat <- addMeta(cellchat, meta = sc_data@meta.data)
cellchat <- setIdent(cellchat, ident.use = "celltype")

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
netVisual_circle(cellchat@net$count, vertex.weight = T)
```

### 4.6 调控因子分析（SCENIC）
```R
# 提取免疫细胞表达矩阵
expr_mat <- as.matrix(GetAssayData(immune_clusters, slot = "data"))
# SCENIC 分析见独立章节（TF 调控网络分析）
```

---

## 5.结果解读与可视化
| 分析结果              | 含义              |
| ------------------ | ------------------------ |
| **UMAP / tSNE**    | 展示免疫细胞亚群分布      |
| **Marker heatmap** | 各免疫亚群特征基因表达     |
| **GSVA 热图**        | 各免疫亚群免疫通路活性     |
| **Checkpoint 分布图** | 免疫抑制基因表达水平      |
| **CellChat 网络图**   | 细胞通讯网络与配体-受体强度  |
| **SCENIC 网络图**     | 关键转录因子调控免疫亚群的模式 |


示例可视化：
```R
DimPlot(immune_clusters, group.by = "SingleR_labels")
DoHeatmap(immune_clusters, features = top10_markers)
netVisual_bubble(cellchat, sources.use = "T_cells", targets.use = "Macrophages")
```
---

## 6.总结
肿瘤免疫微环境分析（TME） 是连接转录组数据与免疫生物学解释的关键环节。通过整合 scRNA-seq、CellChat、GSVA、SCENIC 等模块，可在单细胞层面实现：
- 精细化免疫细胞分型
- 功能状态与信号通路解析
- 免疫逃逸与细胞通讯机制揭示
- 关键调控因子与潜在治疗靶点预测
> 肿瘤免疫微环境分析是通过转录组等组学数据，在细胞层面上系统揭示**肿瘤–免疫系统的动态交互与调控机制**，是当前精准肿瘤免疫研究的核心方向之一。

---

## 7.常见问题解答（FAQ）
**1.TME 分析可以用 bulk RNA-seq 吗？**

可以，但只能估算免疫细胞比例或群体表达特征，无法解析单细胞异质性。常用工具：`CIBERSORTx`、`EPIC`、`xCell`、`MCP-counter`。

**2.scRNA-seq 数据做 TME 分析的优势是什么？**

可以识别免疫细胞亚群（如 CD8+ T 细胞耗竭状态、M1/M2 巨噬细胞）；可分析细胞状态、配体-受体通讯和转录调控网络；可结合空间信息揭示免疫细胞与肿瘤细胞的空间关系。

**3.样本量少能做 TME 分析吗？**

对于 bulk RNA-seq，样本少会影响统计功效和免疫细胞比例估算准确性。对于 scRNA-seq，单个样本细胞数多可以部分弥补样本少的问题，但多样本仍有助于结果稳健性。

**4.TME 分析结果如何验证？**

可以使用实验验证和交叉数据集验证。实验验证：免疫组化（IHC）、流式细胞术（FACS）验证细胞类型和比例。交叉数据验证：使用独立数据集或不同平台验证结果一致性。

**5.免疫细胞比例估算准确吗？**

bulk RNA-seq 方法依赖基因表达参考矩阵，估算精度受限于样本类型、肿瘤异质性及表达噪声。scRNA-seq 可提供单细胞水平分辨率，精度更高。

**6.TME 分析是否需要去除肿瘤低质量样本？**

是的，样本质量差会影响免疫细胞比例估算或亚群识别。QC 指标可包括 <u>RNA 完整性、细胞数量、线粒体比例、样本批次信息</u>等。

**7.可以整合多组学数据进行 TME 分析吗？**

可以，整合突变组学、甲基化、代谢组数据，可深入理解免疫逃逸和肿瘤-免疫调控机制。当然如果分析之后有实验验证将会更加完美。

**8.TME 分析能预测免疫治疗反应吗？**

部分工具或特征（如 CD8+ T 细胞比例、免疫耗竭分数、免疫检查点表达）可作为免疫治疗预测指标，但需结合临床验证。

---