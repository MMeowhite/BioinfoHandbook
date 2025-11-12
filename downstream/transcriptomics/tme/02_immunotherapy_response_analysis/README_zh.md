# 免疫检查点与免疫治疗反应预测（Immunotherapy Response Analysis）
免疫治疗（Immunotherapy）特别是免疫检查点抑制剂（Immune Checkpoint Inhibitors, ICIs）已经成为癌症治疗的重要手段。然而，不是所有患者都能获益，因此预测患者对免疫治疗的反应具有重要临床价值。本章介绍免疫检查点分析与免疫治疗反应预测的基础方法和流程。

---

## 1.什么是免疫检查点与免疫治疗反应预测？
首先需要了解一个重要的概念：免疫检查点。免疫检查点是指**免疫系统中的抑制性分子，可调控 T 细胞活性、防止自身免疫**。常见的免疫检查点：
- PD-1（Programmed Cell Death Protein 1）
- PD-L1（Programmed Death Ligand 1）
- CTLA-4（Cytotoxic T-Lymphocyte Associated Protein 4）
- LAG3、TIM3、TIGIT 等

作用：肿瘤细胞通过上调这些抑制分子逃避免疫攻击。

而免疫治疗反应预测则是**利用基因表达、突变负荷或免疫细胞特征预测患者对免疫检查点抑制剂（如PD-1/PD-L1或CTLA-4抗体）的响应**。常用指标如下：
- 免疫检查点表达水平（PD-1/PD-L1/CTLA-4）
- 肿瘤突变负荷（TMB）
- 微卫星不稳定性（MSI）
- 免疫细胞浸润（TILs）
- 免疫特征评分（Immune Score / TIDE score）

---

## 2.常用工具与数据库
### 分析与可视化工具（R / Python）：
- **CIBERSORTx / xCell / MCP-counter**：免疫细胞比例估算
- **TIDE (Tumor Immune Dysfunction and Exclusion)**：预测免疫治疗反应
- **ImmuneCellAI**：T 细胞亚群与免疫治疗反应预测
- **ESTIMATE**：计算免疫评分和基质评分
- **ggplot2 / pheatmap / ComplexHeatmap**：免疫检查点可视化

### 常用的数据库
- **TCGA**：癌症基因组与表达数据，适合免疫相关分析
- **cBioPortal**：临床关联数据与突变信息
- **MSigDB**：免疫相关基因集（例如 Hallmark Immune Genes）
- **ImmPort / ICGC**：免疫分子与临床响应数据

---

## 3.数据输入与准备
输入数据：
- 基因表达矩阵（bulk RNA-seq 或 scRNA-seq）
- 样本临床信息（是否接受免疫治疗、治疗反应、分期等）
- 免疫细胞比例或评分矩阵（可通过 CIBERSORTx / xCell / MCP-counter 得到）

注意事项：
- bulk RNA-seq 使用原始 counts 或归一化数据均可，但需要保证批次一致性
- scRNA-seq 可以计算特定细胞亚群的免疫检查点表达

---

## 4.分析流程
### 4.1 免疫检查点表达分析
```R
# 选择常用免疫检查点基因
checkpoint_genes <- c("PDCD1","CD274","CTLA4","LAG3","HAVCR2","TIGIT")

# 提取表达矩阵
checkpoint_expr <- expr_matrix[checkpoint_genes, ]

# 可视化热图
library(pheatmap)
pheatmap(checkpoint_expr, scale="row", cluster_rows=TRUE, cluster_cols=TRUE,
         main="Immune Checkpoint Expression")
```

### 4.2 免疫评分与细胞浸润估算
```R
# 使用 ESTIMATE 计算 Immune Score
library(estimate)
# 输入已格式化表达矩阵
immune_scores <- estimateScore(expr_matrix)
# 可视化免疫评分分布
boxplot(immune_scores$ImmuneScore ~ clinical_info$response,
        main="Immune Score vs Response", ylab="Immune Score")
```

### 4.3 TIDE 分析
```R
# 安装 TIDE R 包或在线使用
# 预测免疫治疗反应
# 输出 TIDE score: 高分 -> 低反应概率，低分 -> 高反应概率
```

### 4.4 差异免疫检查点分析
```R
library(limma)
design <- model.matrix(~ response, data = clinical_info)
fit <- lmFit(checkpoint_expr, design)
fit <- eBayes(fit)
topTable(fit)
```
> 可筛选出在响应组与非响应组之间显著差异的免疫检查点基因

### 4.5 可视化
可视化的类型分为多种：热图（checkpoint expression across samples），箱线图（响应 vs 非响应），小提琴图（scRNA-seq各细胞群的免疫检查点表达）。
```R
# 箱线图示例
library(ggplot2)
ggplot(data=checkpoint_expr_long, aes(x=response, y=expression, fill=response)) +
  geom_boxplot() + facet_wrap(~gene) + theme_bw()
```

---

## 5.结果解读与导出
### 5.1 免疫检查点表达
| 分析结果    | 含义                  |
| ------- | ------------------- |
| 免疫检查点热图 | 展示各样本/细胞群的基因表达水平    |
| 差异免疫检查点 | 显示在响应组与非响应组中上调/下调基因 |

### 5.2 免疫评分
| 分析结果                          | 含义                    |
| ----------------------------- | --------------------- |
| Immune Score / ESTIMATE Score | 高分表示免疫浸润丰富，潜在更易响应免疫治疗 |
| TIDE Score                    | 高分表示可能免疫抑制，反应低        |

### 5.3 数据导出
```R
# 导出免疫检查点表达
write.csv(checkpoint_expr, "immune_checkpoint_expression.csv", row.names=TRUE)

# 导出免疫评分与TIDE结果
write.csv(cbind(clinical_info, immune_scores, TIDE_score), "immune_response_scores.csv", row.names=TRUE)
```

## 6.总结
免疫检查点分析与免疫治疗反应预测帮助研究者：
- 识别关键抑制性免疫分子
- 估算肿瘤免疫微环境活性
- 筛选潜在免疫治疗敏感患者
- 提供免疫治疗相关生物标志物候选

> 建议结合 TME 分析、差异基因分析或 scRNA-seq 分析的结果进行综合解读，增加预测准确性。

---

## 7.常见问题解答（FAQ）
**1.可以直接用 TPM/FPKM 做分析吗？**

可以，但需保证批次一致性；对差异分析建议使用归一化计数或 log2(TPM+1)。
**2.样本量少能预测免疫治疗反应吗？**

样本少可能导致统计功效不足，建议多样本或结合多数据集验证。

**3.scRNA-seq 可以做免疫治疗预测吗？**

可以，通过计算免疫检查点在各细胞群表达及免疫评分，但预测模型一般基于 bulk RNA-seq。

**4.TIDE score 越高越好吗？**

不，TIDE score 越高说明免疫逃逸可能性越大，反应概率越低。

**5.是否可以结合突变数据（TMB/MSI）？**

可以，结合 TMB、MSI 与免疫评分可提高预测准确性。

---
