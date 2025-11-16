# 生信流程（bioinformatics pipeline / workflow）
## 1.什么是生信流程？
生信流程工具是用于**自动化、规范化、可重复地**执行生物信息学分析的一类软件或框架。
它们的核心目的是把生信分析从“手动敲命令”变成“自动化、一键运行、可复现的流程”。

### 1.1 生信流程工具的目的
由于生物信息学分析通常包含很多步骤，例如：原始数据质控 -> 读段比对 -> 去重复、过滤低质量片段 -> 计数与定量 -> 差异分析 -> 注释、可视化 -> 下游统计分析

如果这些步骤手动操作，会出现以下的问题：
- **容易出错**：参数写错、文件路径错、版本不一致等。
- **难以重复**：过几个月换电脑、换环境就跑不起来。
- **难以扩展**：手动跑几十个样本还行，几百上千就不可能。
- **难以共享**：别人想要复现你的结果很困难。
- **耗费时间**：很多步骤其实可以并行，但手动很难操作。

所以为了让生信分析更加**自动化、标准化、可复现、可扩展（大规模）**，我们需要 workflow 工具。

### 1.2 生信流程工具到底做了什么？
可以把生信流程工具理解为一个帮你**管理所有分析步骤的调度系统**：
- **各步骤的串联**：例如QC → 比对 → 去重复 → 计数 → 差异分析，你只要写一次规则，它就会自动判断什么先做、什么后做。
- **自动处理依赖关系**：只有比对完成后，才能计算覆盖度。
- **自动并行化**：如果你有 100 个样本，它可以帮你 同时跑多个样本，节省大量时间。
- **错误恢复**：跑到一半服务器崩了不需要从头开始，它会自动恢复到最近的步骤继续跑。
- **环境管理**：很多工具支持 conda、Docker、Singularity，保证在任何机器上都能跑出完全一样的结果。
- **审计记录**：每一步运行了什么命令、用了什么参数，全部自动记录，便于论文提交或结果复现。


> **生信 workflow 工具就是帮你把一堆命令，变成一个能自动跑的生产线**。你只需要告诉它“产线怎么搭”，它就会：自动抓取原料（输入文件）、按顺序运行机器（分析步骤）、检查错误、保存所有过程、最后把成果交给你（输出结果）。

### 1.3 适用人群
适用多种目的的人群，包括但不限于：
- 做 RNA-seq / ATAC-seq / WGS / metagenomics 等项目的研究者
- 生信新手，希望减少错误、提高效率
- 实验室负责人，需要维护长期稳定、可复现流程
- 生信工程师，需要批量处理大量数据
- 云计算、生信平台开发人员

> 适用于需要稳定、可重复、自动化处理生物信息数据的科研人员、生信分析师和实验室团队。

---

## 2.常见生信流程 / 工作流
### 2.1 常见的生信流程 / 工作流工具详细介绍
生信工作流工具主要用于：自动化、可重复、可扩展地运行多步骤的分析流程。下面按目前使用最广、生态最成熟的工具进行介绍。
#### 2.1.1 Nextflow
该工具是目前全球生信领域最主流的工具，它是基于 Groovy DSL 的 workflow 工具，广泛用于大规模、高复杂度、可移植的生信分析。

特点：
- 极强的并行能力：自动调度任务，能高效跑大规模数据
- 容器支持完善：Docker / Singularity / Conda / Podman 全支持
- 跨平台：可在本地、HPC、AWS、Google Cloud 上无缝运行
- 容错性强：失败步骤自动重试，支持 checkpoint
- 生态强大：与 nf-core 配套（高质量标准流程库）

适用场景：
- 大规模测序项目（WGS/WES/ATAC/RNA-seq）
- 实验室长期维护的标准化 pipeline
- 云计算环境
- 需要可复现、高可靠工作流的团队

#### 2.1.2 nf-core（Nextflow 官方社区流程库）
基于 Nextflow 的 标准化生信流程集合，由社区维护。

特点：
- 每个流程有统一规范、严格测试、自动化 CI
- 所有流程均支持容器
- 文档完善
- 大量成熟流程：RNA-seq、ATAC-seq、ChIP-seq、WGS、metagenomics 等

适用场景：
- 想“开箱即用”的用户
- 想快速搭建高质量标准流程的实验室

#### 2.1.3 Snakemake（Python 用户最友好的 workflow）
基于 Python 的 DSL，语法类似 Makefile + Python。

特点：
- 容易上手（尤其对会 Python 的人）
- 支持 conda、mamba、Docker、Singularity
- HPC、云、集群都能跑
- 可视化 DAG、自动追踪依赖
- 社区成熟、教程丰富

优势：
- 代码直观，可读性强
- 很适合科研人员自行开发 pipeline

适用场景：
- 中等规模分析
- Python 用户
- 学术实验室内部流程管理

#### 2.1.4 WDL + Cromwell（Broad Institute 标准）
WDL（Workflow Description Language）是一种结构化语言；Cromwell 是执行它的引擎。

特点：
- 语法清晰，结构化强
- Broad 使用（GATK 官方流程均采用）
- 云兼容性强（尤其是 Google Cloud）
- 支持任务级别的资源定义（CPU/内存）

适用场景：
- 大型基因组中心
- 使用 GATK、DRAGEN 等标准工具链
- 跨机构、标准化流程构建


#### 2.1.5 CWL（Common Workflow Language）
一种通用规范，使用 YAML/JSON 定义工具和 workflow。

特点：
- 强调**跨平台可移植性**
- 独立于运行引擎（Toil、Cromwell、Arvados 等都能跑）
- 严格的规范化，适合长期维护

适用场景：
- 多团队合作
- 需要高度可复现（reproducibility）
- 大项目、公开数据平台（如 GA4GH）

#### 2.1.6 Galaxy（图形界面 Workflow 平台）
Galaxy是一个 Web GUI 平台，不需要编程即可创建 workflow。

特点：
- 非常适合生物学家
- 有上千个内置工具
- 易于分享流程
- 能搭建本地或云端 Galaxy 服务器

适用场景：
- 不会写代码的用户
- 教学、课程、培训
- 小/中规模分析
- 原型验证（prototype）

#### 2.1.7 Airflow / Argo / Luigi（通用 workflow 工具）
虽然不是为生信专门设计，但也常被用于生信平台开发。

特点：
- Python（Airflow/Luigi）或 Kubernetes（Argo）
- 更适合企业级任务调度
- 用于构建内部的数据平台

适用场景：
- 生信 SaaS 平台
- 企业级 LIMS/生信系统
- 多团队共享数据管道

#### 2.1.8 AI 生成工作流工具（新兴方向）
近年来开始出现：
- 自动从命令历史生成 Snakemake（如 Snakemaker）
- LLM自动生成 Nextflow 或 WDL 原型
- 自动优化步骤依赖与资源调度的研究工具

目前仍在快速发展阶段，但未来潜力大。

### 2.2 总结
下面的对比表总结了当前生物信息学领域最常用的工作流（pipeline/workflow）管理工具。不同工具在设计理念、描述语言、执行方式、可扩展性以及适用场景上存在明显差异。本表重点比较了它们的语言/DSL、容器支持、在集群与云环境中的兼容性，以及在实际项目中的优势、限制与典型应用场景。

| 工具                       | 类型             | 语言 / DSL        | 容器支持                         | HPC/集群 | 云平台支持                  | 优势                          | 劣势                    | 典型应用场景                    |
| ------------------------ | -------------- | --------------- | ---------------------------- | ------ | ---------------------- | --------------------------- | --------------------- | ------------------------- |
| **Nextflow**             | 工作流框架          | Groovy DSL      | Docker / Singularity / Conda | ✔      | AWS / GCP / Azure      | 并行能力强；容错；跨平台；社区生态强（nf-core） | 学习曲线略陡；Groovy 不够常见    | 大规模测序项目、核心设施、可复现流程        |
| **nf-core**              | 标准流程库          | Nextflow        | 内置容器                         | ✔      | ✔                      | 高质量标准流程；规范统一；开箱即用           | 灵活度比自己写稍低             | 快速搭建成熟流程；标准化分析            |
| **Snakemake**            | 工作流框架          | Python 风格 DSL   | Conda / Docker / Singularity | ✔      | AWS / GCP / Kubernetes | 上手快；规则简单；可读性好；适合科研          | 大超大规模任务效率稍逊           | 中小型项目；实验室内部分析流程           |
| **WDL + Cromwell**       | 工作流语言 + 执行引擎   | WDL（类 JSON 结构化） | Docker / Singularity         | ✔      | 强（GCP 特别好）             | Broad 官方支持；结构清晰；大型项目常用      | 生态相对 Nextflow 小；语法有约束 | GATK 流程；大型基因组中心           |
| **CWL**                  | 工作流标准          | YAML / JSON     | Docker                       | 依赖执行器  | 多平台                    | 高可移植性；跨团队合作好；标准化程度高         | 写法繁琐；学习成本高            | 标准化、长期维护、多机构合作            |
| **Galaxy**               | 可视化工作流平台       | 无需代码            | 内置                           | 集群可配   | 可搭建云版本                 | GUI 友好；适合初学者；工具丰富           | 不适合特别大规模数据；灵活性较弱      | 教学、无编程背景用户、小型/中型项目        |
| **Airflow**              | 通用调度系统         | Python          | 可整合                          | ✔      | 强                      | 企业级调度；适合平台开发                | 非专用生信；数据流灵活性不足        | 企业生信平台、pipeline 调度        |
| **Argo Workflows**       | Kubernetes 工作流 | YAML            | Docker                       | 基于 K8s | 强                      | 云原生；可扩展性强；适合大规模并行           | 依赖 Kubernetes；成本高     | 云平台生信服务、工业级项目             |
| **Luigi**                | 通用 workflow    | Python          | 可整合                          | ✔      | 有限                     | 简单任务调度方便                    | 功能较弱；不适合复杂生信 workflow | 小型调度、内部工具链整合              |
| **Toil**                 | 执行引擎           | Python API      | Docker                       | ✔      | ✔                      | 支持 WDL/CWL；大规模任务友好          | 配套生态小，入门资料少           | 多格式 workflow 执行、跨云执行      |
| **Snakemaker**（新兴 AI 辅助） | 自动流程生成工具       | Snakemake       | —                            | —      | —                      | 自动把命令历史转成 Snakemake 流程      | 仍在发展中                 | 快速将 ad-hoc 分析转换为 workflow |


> 总体来说，**Nextflow、Snakemake、WDL/Cromwell 和 CWL 是目前生信分析中最主流的四大工作流体系**，各自拥有成熟生态与广泛用户基础；Galaxy 提供了可视化界面，更适合入门和教学场景；而 Airflow、Argo、Luigi 等工具虽然并非为生信设计，但在生信平台建设与企业级数据调度中也具有重要角色。

--- 

## 3. 如何选择这些工具？

在生物信息学分析中，常见的工作流工具可分为**命令行工具（CLI）**、**图形界面平台（GUI）**和**任务调度/平台工具**。下表总结了主要工具、使用方式特点，并附上跳转到对应详细文档的链接，便于深入学习。

### 3.1 工具分类与使用概览

| 工具 | 类型 | 使用方式特点 | 详细文档 |
|------|------|----------------|-----------|
| **Nextflow** | CLI 工作流框架 | 写 Groovy DSL → `nextflow run`；自动并行化、容器化、容错 | [Nextflow](nextflow/README.md) |
| **nf-core** | 标准流程库 | 基于 Nextflow 的标准流程 → `nextflow run <pipeline>`；开箱即用、社区维护 | [nf-core](nf_core/README.md) |
| **Snakemake** | CLI 工作流框架 | 写 Snakefile → `snakemake -j N`；Python 风格 DSL，上手快，依赖自动解析 | [Snakemake](snakemake/README.md) |
| **WDL + Cromwell** | CLI 工作流语言 + 执行器 | 写 WDL + JSON → `cromwell run`；适合标准化、跨平台、云环境 | [WDL + Cromwell](wdl_cromwell/README.md) |
| **CWL** | CLI 工作流标准 | 写 YAML → `cwltool` 运行；强调可移植性、跨平台、长期维护 | [CWL](cwl/README.md) |
| **Galaxy** | GUI 平台 | 上传文件 → 点击运行 → 得到结果；无需编程，适合初学者与教学 | [Galaxy](galaxy/README.md) |
| **Airflow** | 调度系统 | 写 Python DAG → 调度运行；适合企业级任务、平台集成 | [Airflow](airflow/README.md) |
| **Argo Workflows** | Kubernetes 工作流 | 写 YAML → `argo submit`；云原生、高扩展性、适合大规模并行 | [Argo](argo/README.md) |

> **概览总结**：
> - **Nextflow / Snakemake**：适合科研人员 CLI 使用，自动化、可复现、支持容器和 HPC/云环境。
> - **nf-core**：标准化的 Nextflow pipeline，快速搭建成熟分析流程。
> - **WDL + Cromwell / CWL**：强调标准化、跨平台可复现，适合大型项目或多机构协作。
> - **Galaxy**：图形界面平台，适合无编程背景用户、小型分析和教学。
> - **Airflow / Argo**：任务调度系统或云平台工具，适合平台级部署和大规模项目管理。


### 3.2 使用建议

1. **项目规模**：
    - 小型 / 中型科研分析：Snakemake 或 Galaxy 最合适。
   - 大型测序中心 / 云环境：Nextflow + nf-core / WDL + Cromwell。
   - 企业或平台级调度：Airflow / Argo。

2. **团队技能**：
   - 熟悉 Python：Snakemake、Airflow。
   - 熟悉 Groovy / CLI：Nextflow、nf-core。
   - 编程不熟悉：Galaxy。

3. **可复现性**：
   - 需要严格可复现：Nextflow + nf-core、WDL/CWL。
   - 需要快速原型 / 教学：Galaxy。

4. **容器与依赖管理**：
   - Nextflow / Snakemake / WDL / CWL 都支持 Docker 或 Conda。
   - Galaxy 内置工具和依赖，用户无需额外管理。
   - Argo / Airflow 可集成容器，但需配置 Kubernetes 或调度环境。


### 3.3 一句话快速记忆

- **Nextflow**：Groovy DSL → `nextflow run`  
- **nf-core**：标准 Nextflow pipeline → `nextflow run <pipeline>`  
- **Snakemake**：Snakefile → `snakemake -j N`  
- **WDL/Cromwell**：WDL + JSON → `cromwell run`  
- **CWL**：YAML → `cwltool`  
- **Galaxy**：上传 → 点击运行 → 得到结果  
- **Airflow**：Python DAG → 调度运行  
- **Argo**：YAML → `argo submit`

> 总之，常用的生信工作流工具可以根据不同的工具类型、使用方式、个人喜好和适用场景进行选取。
---
