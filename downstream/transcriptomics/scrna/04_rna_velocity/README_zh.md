# RNA Velocity 分析
RNA velocity（RNA 动态分析）是单细胞转录组分析中用于推断细胞状态变化方向与速度的方法。它基于每个细胞中 spliced（成熟转录本） 与 unspliced（未剪接转录本） 的比例，推测转录速率与基因动态变化趋势，从而预测细胞未来的转录状态。
> RNA velocity 可视为轨迹分析（pseudotime）的补充，能揭示细胞发育方向、状态转换趋势及潜在分化路径。

---

## 1.什么是RNA velocity？
核心思想：在转录过程中，基因的 pre-mRNA（unspliced）逐渐被剪接成成熟的 mRNA（spliced）。若某基因在某个细胞中正在“上调”，则该细胞中 unspliced RNA 会暂时升高；若基因“下调”，unspliced RNA 会减少。通过比较 spliced / unspliced 比例，可以估计每个基因的表达变化速度。
分析目标：
- 推断每个细胞的“速度向量”（velocity vector）
- 在UMAP/t-SNE空间中显示细胞发育方向
- 与轨迹分析结合，揭示细胞分化趋势

---

## 2.常用工具
| 工具包                             | 语言       | 功能说明                                               |
| ------------------------------- | -------- | -------------------------------------------------- |
| **velocyto**                    | Python/R | 最早的 RNA velocity 实现，基于 spliced/unspliced 计数矩阵计算速度。 |
| **scVelo**                      | Python   | 目前最常用、功能最强的工具，支持稳态模型与动态模型，可视化强大。                   |
| **Seurat + velocyto.R**         | R        | 与 Seurat 集成良好，可直接在 R 中进行基本可视化。                     |
| **loomR**                       | R        | 用于读取 `.loom` 文件（velocyto 输出格式）。                    |
| **SeuratDisk / SeuratWrappers** | R        | Seurat 与 scVelo 之间格式转换。                            |

---

## 3.数据输入与准备
RNA velocity 需要：
- spliced/unspliced计数矩阵（通常由velocyto从CellRanger输出的BAM文件生成）
- Seurat对象或AnnData对象（用于保存聚类和降维信息）

数据准备方式：
首先使用velocyto run10x从10x输出目录生成.loom 文件，实例代码如下：
```bash
velocyto run10x sample_folder/ refdata-cellranger-GRCh38-3.0.0/
```
然后将`.loom`文件导入R或Python进行后续分析。

---

## 4.分析流程
### 4.1 在 R 中整合 Seurat 与 RNA velocity 数据
```R
library(Seurat)
library(SeuratWrappers)
library(velocyto.R)
library(loomR)

# 读取 Seurat 对象
seu <- readRDS("scRNAseq_SeuratObject.rds")

# 读取 loom 文件
ldat <- ReadVelocity(file = "sample.loom")

# 将 spliced/unspliced 数据添加到 Seurat 对象
seu[["spliced"]] <- CreateAssayObject(counts = ldat$spliced)
seu[["unspliced"]] <- CreateAssayObject(counts = ldat$unspliced)

# 选择 RNA assay 并归一化
DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
```

### 4.2 计算 RNA velocity
```R
library(velocyto.R)

# 计算速度矢量
rvel.cd <- gene.relative.velocity.estimates(
  seu@assays$RNA@counts,
  deltaT = 1,
  kCells = 20,
  fit.quantile = 0.02
)

# 绘制 UMAP 上的 velocity 流向图
show.velocity.on.embedding.cor(
  emb = Embeddings(seu, "umap"),
  vel = rvel.cd,
  n = 200,
  scale = "sqrt",
  cell.colors = ac(seu$seurat_clusters, alpha = 0.5),
  arrow.scale = 2,
  show.grid.flow = TRUE,
  min.grid.cell.mass = 0.5,
  grid.n = 40,
  main = "RNA Velocity on UMAP"
)
```

### 4.3 附：在 Python（scVelo）中分析
```python
import scvelo as sc
import scanpy as sca

# 读取 loom 文件
adata = sc.read_loom("sample.loom")

# 预处理与标准化
scv.pp.filter_and_normalize(adata, min_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# 计算 RNA velocity
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)

# 可视化
scv.pl.velocity_embedding_stream(adata, basis='umap', color='clusters')
```

## 5.结果解读、可视化与导出
### 5.1 核心结果与含义
| 结果对象                         | 含义                          |
| ---------------------------- | --------------------------- |
| **velocity vector**          | 每个细胞在 UMAP 空间的方向与速度         |
| **velocity streamplot**      | 可视化细胞动态流动趋势                 |
| **latent time / pseudotime** | RNA velocity 推测的时间进程        |
| **phase plot**               | 展示 spliced 与 unspliced 表达关系 |

### 5.2 可视化
R中：
```R
# R 中：UMAP 上 velocity 流线
show.velocity.on.embedding.cor(emb = Embeddings(seu, "umap"),
                               vel = rvel.cd,
                               arrow.scale = 3,
                               grid.n = 50,
                               show.grid.flow = TRUE)
```
Python中：
```python
# Python (scVelo)：流线图
scv.pl.velocity_embedding_stream(adata, basis='umap', color='celltype', legend_loc='right')
```
**解读：流线方向表示细胞状态变化方向，速度箭头长度表示变化速度**。
> 可结合轨迹分析结果，验证分化路径一致性。

### 5.3 导出结果
R中：
```R
saveRDS(rvel.cd, file = "RNA_velocity_results.rds")
```
Python中：
```python
adata.write("scvelo_processed.h5ad")
```

---

## 6.总结
RNA velocity 是单细胞分析中用于揭示细胞动态变化的重要方法。
它通过计算基因转录速率（unspliced vs spliced），推断细胞发育方向与未来状态，常用于：
- 验证轨迹分析结果
- 推测分化终点
- 探索细胞状态转换（如激活→分化）
> RNA velocity 结合 Monocle3 轨迹分析，可实现从“静态聚类”到“动态演化”的系统解析。

---

## 7.常见问题解答（FAQ）
**1.RNA velocity 分析需要什么样的数据？**

必须包含 unspliced / spliced 计数信息，通常由 velocyto run10x 生成 `.loom` 文件。

**2.用 R 还是 Python？**

若希望快速结合 Seurat，可用 `velocyto.R`；若要进行动态模型（latent time）、更复杂分析，推荐 scVelo（Python）。

**3.RNA velocity 与轨迹分析有何不同？**

轨迹分析基于细胞间表达相似性（静态）；RNA velocity 基于转录速率变化（动态），可预测未来方向。

**4.为什么 velocity 图中箭头杂乱或不明显？**

原因可能包括：
- 未过滤低质量基因
- 归一化或降维参数不合适
- 高噪声数据（如线粒体高比例样本）

解决方法：调整 `fit.quantile` 或在高变基因上重新拟合。

**5.可以与 Seurat 聚类结果联合可视化吗？**

可以。scVelo 支持将 Seurat 导出的 .h5ad 数据导入，颜色可按 cluster 或细胞类型显示。

**6.RNA velocity 能用于时间序列样本吗？**

可以，但要注意：RNA velocity 推测的是单细胞内部的时间方向，而非采样时间。