# 富集分析

在完成差异基因分析后，常见的下游任务之一是富集分析（Functional Enrichment Analysis），用于将差异基因或基因排序结果映射到已知的生物学注释（如基因本体 GO、通路 KEGG/Reactome、基因集 MSigDB 等），以揭示差异基因的功能倾向和调控网络。本节以与上一节相同的格式、层次和风格，完整且有逻辑地讲解富集分析的原理、工具、数据准备、具体流程、可视化、结果解读、常见问题与注意事项，便于直接整合到你的文档中。

## 1.什么是富集分析

富集分析的目标是判断一组基因（通常是差异基因或某种排序列表）在某些已知功能类别或通路中是否出现得比随机更多（或在整个排序中是否显著富集）。常见类型包括：
- **过度表达富集（Over-Representation Analysis, ORA）**：把显著差异基因作为“兴趣基因集”，与背景基因集合对比，计算每个注释集（如 GO 术语、KEGG 路径）中基因富集的统计显著性（通常基于超几何检验/ Fisher 精确检验）。
- **基因集合富集分析（Gene Set Enrichment Analysis, GSEA）**：使用基因的连续排序（如按 log2FoldChange 或统计量排序），检测整个排序中某个基因集是否系统性地聚集在前端或后端（不依赖差异基因的硬阈值）。
- **富集的变体（Variants）**：如基于权重的GSEA、ssGSEA/GSVA（单样本富集评分）等，用于不同数据与问题场景。

> 富集分析能把千上万的基因列表浓缩为少量可解释的生物学主题（如“免疫反应”、“细胞周期”），是连接基因水平结果和生物学机制的关键桥梁。

---

## 2.常用工具与数据库
### 分析与可视化工具（R）

- **clusterProfiler**：集成 GO、KEGG、GSEA、Reactome 等分析与多种可视化（dotplot、cnetplot、emapplot、ridgeplot 等）。强烈推荐用于基于 R 的工作流。
- **ReactomePA**：Reactome 通路富集与可视化（可与 clusterProfiler 互补）。
- **enrichplot**：clusterProfiler 的可视化扩展（dotplot, emapplot, cnetplot, ridgeplot）。
- **pathview**：将基因表达或富集结果映射到 KEGG 路径图上（图形可视化）。
- **GSVA / ssGSEA（GSVA 包）**：单样本层面的基因集合富集评分，适合样本间比较。
- **fgsea**：高效的 GSEA 实现（基于预计算的无偏置置换算法），适合大基因集或重复运行。
- **DOSE**：疾病富集与可视化（与clusterProfiler联用）。

### 注释数据库与基因集来源

- **GO（Gene Ontology）**：Biological Process / Molecular Function / Cellular Component。
- **KEGG（Kyoto Encyclopedia of Genes and Genomes）**：代谢/信号通路图谱（注意 KEGG 的物种代码）。
- **Reactome**：高质量的通路注释（反应/通路）。
- **MSigDB（Broad）**：Hallmark、C2、C5等预定义基因集合（常用于GSEA）。
- **物种注释**：org.Hs.eg.db、org.Mm.eg.db 等，用于ID转换与注释。

---

## 3.数据输入与准备
### 输入类型（两类主要情形）：
- **1.一组显著差异基因（用于 ORA）**
    - 例如：DESeq2/edgeR 输出经过阈值筛选的上调或下调基因（GENE ID 列表）。
    - 需明确基因 ID 类型（ENSEMBL / ENTREZID / SYMBOL 等）并确保与富集工具/注释数据库一致。
- **2.按某统计量排序的全基因列表（用于 GSEA）**
    - 例如：按stat、log2FoldChange或 -log10(pvalue)*sign(log2FC) 排序的命名向量（names为gene IDs，values为排序统计量）。
    - GSEA要求完整基因排序，不强制二值化。

### 背景基因（gene universe）
背景集合决定了超几何检验的“总体”，应选择与实验可检测到的基因集合一致（例如过滤低表达后用于差异分析的基因集合），而不是全基因组，错误选择会导致假阳性或假阴性风险。

### 基因ID转换
基因的ID分为多种类型，常见的包括Gene symbol，ENSG ID，ENTREZID等，通过R包的函数能够实现基因ID之间的互相转换：Gene symbol <-> ENSG ID<-> ENTREZID。常见的函数：```clusterProfiler::bitr()、AnnotationDbi::select()、biomaRt```等。

> 强烈建议把所有基因映射为ENTREZID（许多R包以ENTREZ为首选），并记录未映射基因数目。

---

## 4.详细流程（以R/clusterProfiler 为例）
### 4.1 环境与数据准备
- **目的**：首先根据差异富集分析得到的表格（包含Gene ID, padj, log2FoldChange等列的表格）进一步筛选，然后讲Gene ID转换为ENTREZID的类型，为后续的富集分析进行相应的数据格式、类型等准备。
```R
# 首先加载相关的R包
library(clusterProfiler)   # 基因富集分析包
library(org.Hs.eg.db)      # 人类注释包，换成对应物种

# 1) 获取显著基因（ORA）
sig_genes <- rownames(subset(res, padj < 0.05 & abs(log2FoldChange) > 1))

# 2) 将 symbol 转为 ENTREZID（示例）
library(clusterProfiler)
gene_map <- bitr(sig_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
entrez_sig <- unique(gene_map$ENTREZID)

# 3) 背景基因（所有用于差异分析的基因）
universe_genes <- rownames(res)
universe_map <- bitr(universe_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
entrez_universe <- unique(universe_map$ENTREZID)
```

- **注意**：

### 4.2 ORA：富集分析
- **目的**：该步骤为正式富集分析的步骤，根据你想要富集的数据库的不同可以选择不同包的不同函数输入即可。

#### 4.2.1 GO富集分析
- **示例代码**：
```R
ego <- enrichGO(gene         = entrez_sig,
                universe     = entrez_universe,
                OrgDb        = org.Hs.eg.db,
                ont          = "BP",            # BP / MF / CC
                pAdjustMethod= "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2,
                readable     = TRUE)            # 返回 SYMBOL 更可读
```
#### 4.2.2 KEGG富集分析
- **示例代码**：
```R
# 若用 ENTREZID
ekegg <- enrichKEGG(gene         = entrez_sig,
                    universe     = entrez_universe,
                    organism     = 'hsa',         # 注意物种代码，如 hsa (human), mmu (mouse)
                    pAdjustMethod= 'BH',
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2)
# 可转换为可读基因名
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
```
#### 4.2.3 Reactome富集分析
- **示例代码**：
```R
library(ReactomePA)
er <- enrichPathway(gene         = entrez_sig,
                    universe     = entrez_universe,
                    organism     = "human",
                    pAdjustMethod= "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2,
                    readable     = TRUE)
```
#### 4.2.4 GSEA富集分析（基于排序向量，例如 log2FC）
- **示例代码**：
```R
# 1. 构建命名向量（ENTREZID 名称）
# 假设 'geneList' 为 data.frame 两列：gene (SYMBOL) 和 stat（例如 log2FC）
geneList <- res$log2FoldChange
names(geneList) <- rownames(res)
# 转换为 ENTREZID 并去 NA（示例）
map_all <- bitr(names(geneList), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# 将重复 SYMBOL 合并（采用最大值或平均，示例取最大）
library(dplyr)
map_df <- as.data.frame(map_all)
gene_df <- data.frame(SYMBOL = names(geneList), stat = geneList)
merged <- merge(map_df, gene_df, by.x="SYMBOL", by.y="SYMBOL")
# 保留最有代表性的值（示例：按 ENTREZID 取最大）
library(tidyr)
library(dplyr)
ranked <- merged %>% group_by(ENTREZID) %>% summarize(stat = max(stat)) 
geneList_entrez <- ranked$stat
names(geneList_entrez) <- ranked$ENTREZID
geneList_entrez <- sort(geneList_entrez, decreasing = TRUE)

# 2. 运行 GSEA（以 GO 为例）
gseaGO <- gseGO(geneList     = geneList_entrez,
                OrgDb        = org.Hs.eg.db,
                ont          = "BP",
                nPerm        = 1000,
                minGSSize    = 10,
                maxGSSize    = 500,
                pAdjustMethod= "BH",
                pvalueCutoff = 0.05,
                verbose      = TRUE)
# 3. 或使用 fgsea（较快）：
library(fgsea)
# 需要 gene sets 列表（例如 MSigDB 或 GO gene sets），然后：
# fgseaRes <- fgsea(pathways = pathwaysList, stats = geneList_entrez, minSize=15, maxSize=500, nperm=10000)
```
> 以上的富集分析为标准的富集分析，如果想要个性定制分析，仍可以通过参考文档进行参数设置进行分析。

### 4.3 可视化（clusterProfiler/enrichplot）
- **目的**：运行完富集分析的计算流程之后，返回的结果需要进一步进行可视化展示。

```R
library(enrichplot)
# dotplot（常用）
dotplot(ego, showCategory=30) + ggtitle("GO BP Enrichment")

# barplot
barplot(ego, showCategory=20)

# cnetplot（基因-通路网络）
cnetplot(ego, showCategory=10, foldChange=geneList_entrez)

# emapplot（富集关系网络）
emapplot(pairwise_termsim(ego))

# ridgeplot（GSEA 结果）
ridgeplot(gseaGO, showCategory=20)

# pathview：可视化 KEGG 路径图（需安装 pathview）
library(pathview)
# 对某条 KEGG 路径（pathway.id）映射表达差异
pathview(gene.data = geneList_entrez, pathway.id = "hsa04110", species = "hsa")
```
> 同样，以上的函数只是为快速展示富集分析可视化提供便捷，如果想要进步调整图样，可以利用ggplot2包自行进行绘图。

### 4.4 导出结果
- **目的**：导出富集分析的表格，一般用系统自带的`write.csv()`函数保存为`.csv`格式即可。

```R
# 保存 ORA 结果表
write.csv(as.data.frame(ego), file="GO_enrichment_ORA.csv", row.names=FALSE)
write.csv(as.data.frame(ekegg), file="KEGG_enrichment_ORA.csv", row.names=FALSE)
# 保存 GSEA 结果
write.csv(as.data.frame(gseaGO), file="GSEA_GO_results.csv", row.names=FALSE)
```

> 当然，也可以利用`openxlsx::write.xlsx()`函数保存为`.xlsx`格式，便于查看。

--- 

## 5.结果解读与下游分析建议
### 5.1 如何阅读富集表
- ID / Description：通路或术语名称。
- GeneRatio / BgRatio：前者表示输入基因集中属于该术语的比例（k / M），后者表示背景中属于该术语的比例（K / N）。
- pvalue / p.adjust / qvalue：统计显著性（多重检验校正后 p.adjust 更可信）。
- geneID：该术语内实际落在输入基因集中的基因（通常以分隔符列出）。
- Count：输入基因集中属于该术语的基因数。

> 优先关注**p.adjust**或**qvalue（FDR 校正）**，并结合GeneRatio与 Count判断生物学意义（例如高度显著但只包含1-2个基因的术语需谨慎）。

### 5.2 从统计显著到生物学解释
- 检查重复或相关术语（常见于 GO），可以利用 pairwise_termsim() 与 emapplot() 合并相关主题。
- 可结合上/下调基因集分别富集分析（上调/下调分开做 ORA 或 GSEA），以揭示方向性生物学过程。
- 用pathview或网络图（Cytoscape）可视化显著通路的基因位置与表达趋势，便于直观解释信号传导或代谢变化。

### 5.3 下游延展分析
- 模块/网络分析：结合 WGCNA 或 PPI（STRING/Cytoscape）识别关键调控模块或 hub 基因。
- 整合多组学：将转录组富集结果与蛋白组、代谢组或表观组数据整合，以强化路径级别的证据。
- 单样本评分：使用 GSVA/ssGSEA 对每个样本打分，进行样本分群或与表型关联分析（如生存分析）。

---

## 6.总结
富集分析把差异基因从基因层级提升到通路/功能层级，是非常重要的生物学解释手段。
- 根据输入类型选择方法：显著基因集合 → ORA；完整排序 → GSEA。两者结合可互补。
- 关键步骤包括：基因 ID 规范化、合适的背景选择、正确的注释数据库、严格的多重检验校正。
- 使用 clusterProfiler / fgsea / ReactomePA / pathview 等工具可以实现从分析到可视化的一体化工作流。
- 结果解读需统计显著性 + 生物学合理性 + 多数据/数据库互证。最后应以实验验证（如 qPCR、功能实验）来确认候选通路与基因的实际生物学作用。

---

## 7.常见问题解答（FAQ）
**1.应该用 ORA 还是 GSEA？**

主要取决于你的分析目的：
- ORA（过度富集）：当你有明确的显著基因列表（阈值筛选后）时使用，简单且解释直观；但它对阈值很敏感（可能遗漏弱但一致的信号）。
- GSEA：使用完整排序信息，能检测到那些在排序中一致偏向前端/后端但未达到硬阈值的基因集合，适合信号较弱或复杂背景的情况。两者常结合使用以互为验证。

**2.背景基因（universe）如何选择？**

背景应为实验中“可检测到”的基因集，例如所有表达量通过过滤的基因或参与差异分析的基因；不要随意使用全基因组（除非实验确实有检测覆盖全基因组），错误的背景会使富集统计失真。

**3.ID 类型不匹配怎么办？**

需要统一ID（建议 ENTREZID）并记录丢失或转换失败的基因数目。使用 clusterProfiler::bitr()、biomaRt 或 AnnotationDbi::select() 做批量转换，若大量基因无法映射，需检查输入类型（ENSEMBL vs SYMBOL）或注释包是否匹配物种。

**4.p值/padj的阈值如何设定？**

常用阈值：`p.adjust < 0.05（BH 校正）`；GSEA 常报告`FDR < 0.25（原 GSEA 文献建议，较宽松）`，但具体阈值应结合生物学背景与后续验证可行性调整。

**5.KEGG中organism 参数如何设置？**

organism为三字母代码（如 `hsa：human`，`mmu: mouse`）。注意 KEGG 数据库权限或接口变化可能影响直接访问，clusterProfiler与pathview提供常用支持；若出错，可使用本地或其他来源的基因集。

**6.GO 结果过于冗余怎么办？**

使用`simplify()（clusterProfiler）`基于相似度去掉冗余GO项；或用 pairwise_termsim() + emapplot() 聚类显示主题；人工合并或聚焦在生物学上更关键的大类也常用。

**7.如何判断富集结果是否可信？**

可以从多个层面进行查看
- 统计上看：多重检验校正后显著（p.adjust），且包含足够数量的基因（Count）。
- 生物学上看：是否与已有文献/已知机制一致；是否在上下游不同分析中复现（例如上调基因集和 GSEA 指向相同通路）。
- 技术上看：背景选择、ID 转换与注释是否无误、样本质控是否充分。

---