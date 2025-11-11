# 单细胞 RNA-seq（scRNA-seq）基础分析

单细胞 RNA-seq 可以解析组织中每个细胞的转录组信息，识别细胞类型、亚群、发育轨迹和状态差异。与 bulk RNA-seq 不同，scRNA-seq 数据通常稀疏且噪声较大，需要专门的方法和工具进行处理和分析。本基础分析主要涵盖核心通用流程，为后续轨迹分析、细胞通讯分析等进阶分析打基础。

## 1. 什么是 scRNA-seq 基础分析？

scRNA-seq 基础分析的目标：
- 数据质控，去除低质量细胞和异常基因  
- 数据归一化与高变基因筛选  
- 降维与可视化（PCA / t-SNE / UMAP）  
- 聚类分析，识别细胞群体  
- 差异表达分析，发现各群体特征基因  
- 初步细胞类型注释  

核心思路：
- 每个细胞的转录组被视为独立样本，基因表达矩阵稀疏  
- 使用归一化和高变基因挑选减小噪声  
- 降维方法提取主要表达模式，用于可视化和聚类  
- 聚类结果结合 marker gene 注释确定细胞类型  

---

## 2. 常用工具与数据库
### 分析与可视化工具（R）

- **Seurat**：最常用的单细胞分析包，支持全流程分析：质控、归一化、高变基因、降维、聚类、差异分析。  
- **SingleCellExperiment**：统一数据结构，兼容 Bioconductor 包。  
- **scran / scater**：归一化、质量控制、检测高变基因。  
- **DoubletFinder / scDblFinder**：检测双细胞。  
- **ggplot2 / pheatmap / ComplexHeatmap / patchwork**：多种可视化方法。  
- **monocle3 / slingshot**：发育轨迹和伪时间分析（进阶模块）。  
- **CellChat / CellPhoneDB**：细胞通讯分析（进阶模块）。  

### 注释数据库与基因集来源
- **CellMarker / PanglaoDB**：已知细胞类型标记基因库  
- **GO / KEGG / Reactome / MSigDB**：功能富集分析与通路注释  
- **Human Cell Atlas (HCA)**：人类参考单细胞表达数据  
- **Ensembl / org.Hs.eg.db**：基因注释与 ID 转换 

---

## 3. 数据输入与准备

输入数据：
- **表达矩阵（counts）**：行基因（gene），列细胞（cell）  
- **样本元信息（metadata）**：每个细胞的来源、分组、批次信息等  

> scRNA-seq 数据稀疏，建议使用原始 counts 进行归一化，不建议直接使用 TPM 或 FPKM。

---

## 4. 分析流程（以 Seurat 为例）

### 4.1 构建Seurat对象
```R
library(Seurat)

# 创建 Seurat 对象
sc_data <- CreateSeuratObject(counts = counts_matrix,
                              project = "scRNAseq",
                              min.cells = 3, min.features = 200)
```

### 4.2 质量控制（QC）

```R
# 计算线粒体基因比例
sc_data[["percent.mt"]] <- PercentageFeatureSet(sc_data, pattern = "^MT-")

# 可视化 QC 指标
VlnPlot(sc_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 过滤低质量细胞
sc_data <- subset(sc_data, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
```

### 4.3 数据归一化与高变基因选择
```R
# LogNormalize
sc_data <- NormalizeData(sc_data, normalization.method = "LogNormalize", scale.factor = 10000)

# 高变基因选择
sc_data <- FindVariableFeatures(sc_data, selection.method = "vst", nfeatures = 2000)
```

### 4.4 数据缩放与降维
```R
# 全基因数据缩放
all.genes <- rownames(sc_data)
sc_data <- ScaleData(sc_data, features = all.genes)

# PCA 降维
sc_data <- RunPCA(sc_data, features = VariableFeatures(object = sc_data))
ElbowPlot(sc_data)  # 选择主成分数量

# UMAP / t-SNE 可视化
sc_data <- RunUMAP(sc_data, dims = 1:20)
DimPlot(sc_data, reduction = "umap", group.by = "sample")
```

### 4.5 聚类分析
```R
sc_data <- FindNeighbors(sc_data, dims = 1:20)
sc_data <- FindClusters(sc_data, resolution = 0.5)
DimPlot(sc_data, reduction = "umap", label = TRUE)
```

### 4.6 差异表达分析
```R
# 查找各 cluster marker genes
cluster_markers <- FindAllMarkers(sc_data,
                                  only.pos = TRUE,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25)
head(cluster_markers)
```
### 4.7 细胞类型注释
- **手动注释**：根据已知 marker gene 对比（CellMarker / PanglaoDB）
- **自动注释**：SingleR / scCATCH / CellAssign

---

## 5. 结果解读、可视化与导出

在完成 scRNA-seq 基础分析后，我们主要得到并保存以下结果：

### 5.1 聚类与降维结果

| 分析结果         | 含义 |
|-----------------|------|
| **UMAP / t-SNE 图** | 低维可视化细胞间相似性，聚类结果直观呈现细胞群体结构 |
| **PCA 图**        | 展示主要变异方向，确定降维使用的主成分数量 |
| **聚类标签（cluster）** | 每个细胞所属的群体，用于后续 marker gene 分析 |

```R
# 可视化 UMAP 聚类
DimPlot(sc_data, reduction = "umap", label = TRUE, pt.size = 0.5)

# PCA 可视化
DimHeatmap(sc_data, dims = 1:10, cells = 500, balanced = TRUE)
```
### 5.2 Marker gene 分析结果
| 分析结果                | 含义                              |
| ------------------- | ------------------------------- |
| **cluster_markers** | 每个聚类的显著上调基因，反映群体特征              |
| **表达热图**            | 展示 marker gene 在各 cluster 的表达分布 |
```R
# Top marker genes 热图
top_markers <- cluster_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(sc_data, features = top_markers$gene) + NoLegend()
```

### 5.3 细胞类型注释结果
| 分析结果                 | 含义                                 |
| -------------------- | ---------------------------------- |
| **cell type labels** | 根据 marker gene 或自动注释得出的细胞类型        |
| **分布图**              | UMAP 或 t-SNE 上显示各细胞类型空间分布，直观展示组织结构 |

### 5.4 数据导出

导出Seurat对象用于后续分析：
```R
saveRDS(sc_data, file = "scRNAseq_SeuratObject.rds")
```

导出marker gene表格：
```R
write.csv(cluster_markers, "scRNAseq_cluster_markers.csv", row.names = FALSE)
```

导出聚类信息与细胞类型注释：
```R
metadata <- sc_data@meta.data
write.csv(metadata, "scRNAseq_cell_metadata.csv", row.names = TRUE)
```
---

## 6. 总结
scRNA-seq基础分析流程逻辑：数据读取 -> 构建Seurat对象 -> 质控与异常值剔除 -> 数据归一化、高变基因筛选 -> 降维（PCA）与可视化（UMAP / t-SNE） -> 聚类识别细胞群体 -> 差异表达分析，marker gene挖掘 -> 初步细胞类型注释
> 该基础分析为后续轨迹分析、细胞通讯分析和功能富集分析提供核心输入。
---

## 7.常见问题解答（FAQ）
**1.可以以直接用 TPM / FPKM 吗？**
不推荐。scRNA-seq 通常使用原始counts进行归一化和分析，TPM/FPKM 会破坏稀疏性。

**2.高变基因数量如何选择？**
一般 1000–3000 个，高变基因用于降维和聚类，数量过少会丢信息，过多增加噪声。

**3.聚类分辨率如何选择？**
`FindClusters`的resolution控制聚类颗粒度，分辨率低->大 cluster，分辨率高->小 cluster，可结合marker gene自行调整。

**4.UMAP/t-SNE 聚类不明显怎么办？**
- 可能原因：QC 未严格、batch effect 明显、或高变基因选择不理想
- 解决方法：SCTransform 或使用批次校正方法（Harmony / Seurat IntegrateData）

**5.样本量较小可以做 scRNA-seq 吗？**

单细胞数量多可以弥补样本少，但建议尽量多样本以提高结果可靠性

**6.可以进一步做功能分析吗？**

可以，marker genes 或 cluster 内基因可进行 GO/KEGG/Reactome 富集分析，也可用于轨迹分析或细胞通讯分析。

---
