# 单细胞 RNA-seq（scRNA-seq）细胞-细胞通讯分析

细胞-细胞通讯（Cell–Cell Communication, CCC）分析旨在揭示细胞间通过**配体–受体相互作用（ligand–receptor interactions）**进行的信号传递网络，从而帮助理解组织中不同细胞群体之间的功能关系与调控机制。

在 scRNA-seq 数据中，我们可以利用细胞类型特异的基因表达特征推断可能的信号通路、配体–受体配对关系，以及关键调控网络。

---

## 1.分析目标

- 推断不同细胞类型间可能存在的配体–受体相互作用  
- 识别信号通路级别的通讯网络（如 Notch、WNT、TNF、TGF-β 等）  
- 分析发出信号的“sender cells”和接收信号的“receiver cells”  
- 可视化细胞通讯强度、网络结构和通路特异性关系  

---

## 2.常用工具

| 工具 | 编程语言 | 特点 |
|------|-----------|------|
| **CellChat** | R | 最常用、支持可视化网络、通路分析与细胞角色分析 |
| **CellPhoneDB** | Python | 基于已知配体–受体数据库，结果稳健，适合多样本比较 |
| **NicheNet** | R | 结合下游靶基因分析，能预测信号影响的转录响应 |
| **iTALK / SingleCellSignalR / Connectome** | R | 轻量化可视化分析 |

> 推荐使用 **CellChat**，功能全面、与 Seurat 高度兼容，支持人/鼠数据库。

---

## 3.数据准备

输入数据：
- Seurat 对象（已完成细胞类型注释）  
- 每个细胞的表达矩阵（log-normalized）  
- 细胞类型（或 cluster）信息  

```R
# 假设已完成基础分析，并含有 cell_type 注释
table(sc_data$cell_type)
```

## 4.分析流程
### 4.1 加载包与准备数据
```R
library(CellChat)
library(patchwork)
library(Seurat)

# 从 Seurat 对象构建 CellChat 对象
data.input <- GetAssayData(sc_data, assay = "RNA", slot = "data")
meta <- data.frame(labels = sc_data$cell_type, row.names = colnames(sc_data))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

# 设置数据库（人或鼠）
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
```

### 4.2 识别配体–受体相互作用
```R
# 预处理
cellchat <- subsetData(cellchat)  # 提取配体–受体基因
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# 构建通讯网络
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
```

### 4.3 计算通讯通路与强度
```R
# 汇总通路信息
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
```

### 4.4 可视化通讯网络
```R
# 细胞群体间通讯网络热图
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE)

# 信号强度（加权网络）
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE)
```

### 4.5 通路特异性分析
```R
# 查看特定信号通路（如 WNT）
netVisual_aggregate(cellchat, signaling = "WNT", layout = "circle")

# 查看不同信号方向
netVisual_heatmap(cellchat, signaling = "TGFb")

# 可视化受体–配体气泡图
pairLR <- extractEnrichedLR(cellchat)
netVisual_bubble(cellchat, sources.use = "Macrophage", 
                 targets.use = "Tcell", signaling = c("TNF","IL2"))
```

### 4.6 分析细胞通讯角色
```R
cellchat <- netAnalysis_computeCentrality(cellchat)

# 识别主要信号发送者与接收者
netAnalysis_signalingRole_network(cellchat, signaling = "TGFb")

# 总体角色气泡图
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
```

---

## 5.结果解读、可视化与导出

| 结果类型                             | 含义              |
| -------------------------------- | --------------- |
| **通讯网络（net.count / net.weight）** | 细胞间交互数量和强度      |
| **信号通路特异网络**                     | 某通路下的细胞–细胞交互结构  |
| **配体–受体对（pairLR）**               | 具体相互作用分子对       |
| **发送者/接收者角色分析**                  | 哪类细胞是主要的信号源或信号靶 |

导出结果：
```R
# 导出通讯网络表
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "cellchat_communication.csv", row.names = FALSE)

# 导出配体–受体信息
write.csv(pairLR, "cellchat_ligand_receptor_pairs.csv", row.names = FALSE)

# 保存 CellChat 对象
saveRDS(cellchat, file = "CellChat_analysis.rds")
```

## 6.总结
细胞-细胞通讯分析帮助理解组织微环境中不同细胞之间的相互调控关系。
**典型分析流程：数据准备 → 构建 CellChat 对象 → 识别配体-受体 → 构建通讯网络 → 分析信号通路与细胞角色 → 可视化与导出结果**。当然也可结合轨迹分析结果，从时间或状态角度探讨通讯变化趋势。

## 7.常见问题解答（FAQ）
**1.CellChat 分析的输入是什么？**

标准化后的表达矩阵 + 细胞类型注释信息。推荐使用Seurat对象直接转换，方便快捷。

**2.可以分析多个样本或条件比较吗？**

可以。CellChat 支持 mergeCellChat() 合并不同条件，并可比较通讯差异。

**3.配体–受体数据库可以自定义吗？**

可以。用户可加载自定义数据库（例如肿瘤特异信号），通过`subsetData()`进行替换。

**4.CellPhoneDB 与 CellChat 有何区别？**

CellChat 功能更全、可视化更强；CellPhoneDB 精度高、适合批量比较分析。

**5.通讯强度与数量有何区别？**

数量（count）反映细胞间通讯事件数量，强度（weight）考虑基因表达水平与交互活性。