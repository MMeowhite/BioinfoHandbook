# 单细胞 RNA-seq（scRNA-seq）轨迹分析

轨迹分析（Trajectory / Pseudotime Analysis）用于研究细胞发育、分化或状态转变的动态过程。它基于单细胞 RNA-seq 数据，利用细胞之间的转录组差异推断潜在的时间顺序或分化路径，从而揭示连续性生物学过程，如干细胞分化、免疫细胞激活或肿瘤进化。

---

## 1.什么是轨迹分析？

轨迹分析的目标：
- 推断细胞在连续生物过程中的发展顺序（pseudotime）  
- 识别关键状态转变和分化分支  
- 找出驱动细胞状态变化的基因（branch-dependent genes）  

核心思路：
- 通过降维和聚类识别细胞群体  
- 根据基因表达相似性构建细胞发展图（graph / tree）  
- 使用算法拟合伪时间路径，映射细胞状态变化  
- 挖掘关键基因或调控因子，分析分支决策点  

---

## 2.常用工具
### 分析工具（R）

- **Monocle3**：常用轨迹分析工具，支持无监督聚类、PCA/UMAP 降维、轨迹推断、分支分析。  
- **Slingshot**：基于聚类结果构建分支曲线，简单高效。  
- **TradeSeq**：基于 Slingshot/Monocle 的统计分析，识别分支特异性基因。  
- **dynverse / dynplot**：整合多种轨迹算法，提供统一接口和可视化方法。  

---

## 3.数据输入与准备

轨迹分析通常基于：
- **标准化或归一化后的表达矩阵**（如Seurat对象、log-normalized counts）  
- **聚类信息**（cluster labels）或降维结果（PCA/UMAP/t-SNE）  

> **建议只保留高质量细胞和高变基因，减少噪声，提高轨迹推断精度。**

---

## 4.分析流程（以 Monocle3 为例）

### 4.1 安装与加载包

```R
if(!require(monocle3)) BiocManager::install("monocle3")
library(monocle3)
```

### 4.2 构建CellDataSet对象
```R
# 从 Seurat 对象导入 Monocle3
cds <- as.cell_data_set(sc_data)

# 设置聚类信息
cds@clusters$UMAP$clusters <- sc_data$seurat_clusters
cds@colData$cell_type <- sc_data$cell_type  # 可选
```

### 4.3 降维与聚类（可选）
```R
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds, resolution = 1e-3)
plot_cells(cds, color_cells_by = "cluster")
```

### 4.4 学习轨迹图与伪时间
```R
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE)

# 指定起始点（根细胞）进行伪时间计算
cds <- order_cells(cds, root_cells = which(cds@colData$cell_type == "Stem"))
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE)
```

### 4.5 分支特异性基因识别
```R
# 根据 pseudotime 或分支分析动态表达基因
deg_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
deg_genes <- deg_genes[deg_genes$q_value < 0.05, ]
head(deg_genes)

# 可视化动态表达基因热图
plot_genes_in_pseudotime(cds[c("GeneA","GeneB","GeneC"),],
                         color_cells_by = "cell_type")
```

---

## 5.结果解读、可视化与导出
### 5.1 伪时间路径与分支

| 分析结果                  | 含义                   |
| --------------------- | -------------------- |
| **轨迹图**               | 展示细胞在低维空间中的发展路径和分支结构 |
| **伪时间（pseudotime）**   | 每个细胞在轨迹中的位置，反映潜在发展顺序 |
| **分支点（branch point）** | 关键状态转换节点，可分析分支特异性基因  |
导出：
```R
# 导出伪时间
pseudotime <- pseudotime(cds)
write.csv(pseudotime, "scRNAseq_pseudotime.csv")
```
### 5.2 动态基因与分支特异性基因
| 分析结果        | 含义                    |
| ----------- | --------------------- |
| **动态基因**    | 在轨迹中随时间显著变化的基因        |
| **分支特异性基因** | 在特定分支上高表达的基因，可能调控分支决策 |
导出：
```R
write.csv(deg_genes, "scRNAseq_trajectory_genes.csv")
```

### 5.3 可视化示例
- UMAP轨迹图，按 pseudotime 上色
- 热图展示动态基因随 pseudotime 变化
- 分支特异基因折线图或点图
```R
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster=FALSE)
plot_genes_in_pseudotime(cds[head(rownames(deg_genes), 20),])
```

---

## 6.总结
轨迹分析帮助研究者：
- 理解单细胞异质性中的连续变化
- 推断细胞发育或状态转变顺序
- 识别关键驱动基因和分支调控因子
> **核心逻辑：降维 -> 聚类 -> 学习轨迹图 -> 伪时间排序 -> 动态基因识别 -> 分支特异性分析**。结果可用于功能富集、调控因子分析或与细胞通讯分析结合。

---

## 7.常见问题解答（FAQ）
**1.scRNA-seq 样本量太少可以做轨迹分析吗？**

样本少但单细胞数量足够可行，伪时间结果可能不稳定。建议多样本、多细胞增加可靠性。

**2.如何选择根细胞（root cells）？**

一般根据生物学知识选择未分化/干细胞群体作为根细胞。若不确定，可使用 Monocle3 自动选择起点。

**3.分支特异性基因如何筛选？**

使用`graph_test`或`differentialGeneTest`等方法，根据q-value或 FDR 筛选显著基因。

**4.Monocle3 与 Slingshot 区别？**

Monocle3 更适合大规模数据，提供完整可视化流程；Slingshot 简单高效，依赖已有聚类结果，便于快速分析。

**5.可以结合功能富集分析吗？**

可以，将动态基因或分支特异基因用于GO/KEGG/Reactome等富集分析，挖掘潜在调控机制。

---