# 天然产物—靶点时序检索的来源感知日期精度政策审计

**作者：** Chenyu Hu  
**单位与联系方式：** 作者需在期刊投稿前确认

## 项目目的

本仓库包含“天然产物—靶点时序检索的来源感知日期精度政策审计”1.0.0
版可重复性软件。研究在正例—未标注语义下，对四个固定检索基线进行事后审计；
它不是一个新预测模型。当前分析入口仅由
`manifests/CURRENT_EXECUTION_CHAIN.json` 指定。

作者侧洁净环境复现已**完成**：19 阶段来源到冻结输入链的 16/16 个输入均与
锁定 SHA-256 一致，12 步分析链通过，7 份修订/v4 验证回执通过，稿件汇总
数字检查为 17/17。独立第三方复现**尚未执行**。证据见
`reproduction/clean_environment/`。

## 数据来源与 NPASS 边界

研究使用 NPASS 2.0/3.0、ChEMBL 31、PubMed 日期元数据、UniProt 映射与
序列、RDKit 2026.03.4 和 MMseqs2 18-8cc5c。来源链接、记录的获取日期、
大小与哈希见 `manifests/source_download_manifest.tsv`。

NPASS 状态为 `LINK_ONLY_NO_REDISTRIBUTION`。仓库不包含 NPASS 原始文件、
标识符子集、结构、活性记录、靶点—化合物—PMID 对应关系，或其他逐行派生物。
用户应从官方页面获取源文件：https://bidd.group/NPASS/downloadnpass.html 。
对应版本论文见 `docs/NPASS_CITATION_AND_PROVENANCE.md`。

This repository is not affiliated with or endorsed by the NPASS database or
its maintainers.

## 环境与运行

记录的 Windows 环境为 Python 3.11.15、NumPy 2.4.6、scikit-learn 1.9.0、
RDKit 2026.03.4、pandas 3.0.3、Biopython 1.87 和 Matplotlib 3.11.0。

```text
conda create --name npass_temporal_release --file environment/conda-win-64-explicit.txt
conda activate npass_temporal_release
```

较轻量的跨平台环境：

```text
conda env create -f environment/environment.yml
conda activate npass_temporal_release
```

MMseqs2 不随仓库分发，须从官方项目安装 18-8cc5c 版。

仅检查官方端点：

```text
python scripts/download_sources.py --manifest manifests/source_download_manifest.tsv --check-only
```

把哈希锁定的源文件下载到被忽略的本地目录：

```text
python scripts/download_sources.py --manifest manifests/source_download_manifest.tsv --download --download-dir work/sources --allow-large
```

阅读 `docs/SOURCE_AT_ACQUISITION.md`、`docs/NPASS_LINK_ONLY_POLICY.md` 和
`reproduction/RECONSTRUCTION_REPORT.md`。将路径模板另存为不会提交的
`configs/reproduction_paths.local.json`，填入本地路径后运行：

```text
python scripts/rebuild_analysis.py --config configs/reproduction_paths.local.json
```

该命令对 16 个本地输入执行 fail-closed 哈希门禁。可在一个全新的目录中建立
历史执行布局（不会复制任何受限输入）：

```text
python scripts/materialize_execution_layout.py --output work/execution-layout
```

随后按 `reproduction/protocol/execution-runbook.md` 和
`reproduction/SOURCE_TO_FROZEN_PROVENANCE.json` 执行。缺少前置文件或哈希
不符时，流程会立即停止。PubMed 动态请求会记录 ID
清单哈希、UTC 请求时间、批次、重试及响应 receipt；gzip 采用固定元数据和
稳定排序。

无需第三方真实记录的安全检查：

```text
python scripts/run_reproduction.py --mode verify-chain
python scripts/run_reproduction.py --mode smoke
python -m unittest discover -s tests -v
```

## 预期结果与限制

冻结评价包含 4,123 个候选人类单蛋白靶点、222 个查询、358 个后记录严格
A/B 关系和 4,990 个历史关系。保守的仅精确到日政策纳入 20,647 行中的
13,885 行；区间删失将 20,455 个日期已解析行全部判为确定早于截止日，并
改变 141 个证据等级，但未改变历史成员或冷启动分母。精确数字一致性、图和
表重建以发布审计 receipt 为准。

“后记录”不等同于首次生物学发现；未记录配对是未标注而非负例。审计为事后、
结果已知、作者运行且非独立。PubMed 与 UniProt 服务会动态变化，必须使用
receipt 判断是否需要恢复归档响应。ChEMBL 31 约需 4.51 GB 压缩空间及
23.75 GB 解压空间。尚未发生时，不声称完成独立第三方复现。

## 许可证边界

项目原创软件采用 MIT，Copyright (c) 2026 Chenyu Hu。只有在数据集 manifest
中明确通过审查的非重建性汇总文件采用 CC BY 4.0。上述许可证均不覆盖 NPASS、
ChEMBL、PubMed、UniProt 及其记录、第三方二进制、软件包、权重或缓存。详见
`LICENSE_SCOPE.md` 和 `THIRD_PARTY_NOTICES.md`。

在公开对象实际存在并完成验证前，本文件不填写仓库 URL、commit、release URL
或 DOI。
