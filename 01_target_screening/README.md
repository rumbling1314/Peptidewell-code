# 32 菌种特异表面蛋白筛选 —— 脚本说明

本目录包含从**参考基因组下载**到**种间特异性验证**的完整脚本链（4 个脚本 + 本文档）。

全部脚本为纯 Python 3 / Bash 实现，**不依赖第三方 Python 包**（仅需系统命令行工具），可在 Linux / WSL / macOS 下运行。

---

## 一、环境依赖

| 依赖 | 版本要求 | 用途 | 安装 |
|---|---|---|---|
| **NCBI Datasets CLI** | v14+ | 下载 RefSeq 参考基因组蛋白序列 | https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/ |
| **BLAST+** | 2.10+ | 本地蛋白序列比对（`makeblastdb` / `blastp`） | `conda install -c bioconda blast` 或 https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ |
| **Python** | 3.6+ | 运行脚本（标准库即可） | 系统自带 |
| `unzip` / `curl` | - | 解压与网络请求 | 系统自带 |

可选：NCBI API Key（用于第 3 步在线 BLAST 提频，非必需）
- 申请：https://www.ncbi.nlm.nih.gov/account/ → Account Settings → API Key Management
- 提供方式（任选其一）：
  - 环境变量：`export NCBI_API_KEY=你的key`
  - 文件：`echo '你的key' > ~/.ncbi_api_key && chmod 600 ~/.ncbi_api_key`

---

## 二、脚本流程

```
① 01_download_proteomes.sh     下载 32 菌蛋白组（RefSeq，每菌一个 .faa）
        ↓
② 02_blast_specificity_v2.py   32 菌两两互比（每菌 vs 其余 31 菌合并库）
        ↓
③ 04_rank_specificity_graded.py 按种间一致性分级（S/A/B/C 级）
        ↓
④ 03_batch_blast_ncbi.py       候选蛋白提交 NCBI nr 库终验（在线 BLASTp）
```

---

## 三、各脚本用法

### ① 01_download_proteomes.sh —— 下载 32 菌蛋白组

```bash
bash 01_download_proteomes.sh proteomes
```

- 输入：无（脚本内置 32 菌株的菌名 + RefSeq 组装编号 GCF_xxx）
- 输出：`proteomes/<菌名>.faa`（每菌一个蛋白 FASTA，如 `bacillus_subtilis.faa`）
- 特性：**断点续跑**（已存在的 .faa 自动跳过）
- 说明：使用 `datasets download genome accession GCF_xxx --include protein` 精确下载目标菌株（保证可复现）

### ② 02_blast_specificity_v2.py —— 32 菌互比

```bash
python3 02_blast_specificity_v2.py proteomes blast_output_32strains
```

- 参数：`<蛋白组目录> <输出目录>`
- 逻辑：对每个目标菌（32 轮）：
  1. 将其余 31 菌的蛋白组**合并**为一个 BLAST 库（`makeblastdb -dbtype prot`）
  2. 用目标菌蛋白组查询该合并库：
     ```
     blastp -query <目标菌>.faa -db <31菌合并库> \
            -evalue 1e-3 -outfmt 6 -max_target_seqs 5 -num_threads 4
     ```
  3. 对每个蛋白取**最佳匹配**（最高 bitscore），输出排序表
- 输出：`<输出目录>/<菌名>_specificity_ranked.tsv`
  ```
  protein_id  identity_%  evalue  bitscore  matched_strain  subject_id  aln_length
  ```
- 特性：**断点续跑**（已有结果的菌自动跳过）；自动清理中间文件

### ③ 04_rank_specificity_graded.py —— 分级（可选，筛选用）

```bash
python3 04_rank_specificity_graded.py blast_output_32strains proteomes graded_specificity
```

- 参数：`<互比结果目录> <蛋白组目录> <输出目录>`
- 逻辑：按每个蛋白与最近缘菌种的一致性分级
  - **S 级**：无任何 hit（种内特有）
  - **A 级**：identity < 30%（高度特异）
  - **B 级**：identity 30–40%（中等特异，需后续差异位点分析）
  - **C 级**：identity > 40%（同源度过高，默认淘汰）
- 输出：`<输出目录>/<菌名>_graded.tsv` + `<菌名>_SA_grade.faa`（S+A 级蛋白序列，供下游预测使用）

### ④ 03_batch_blast_ncbi.py —— NCBI nr 库在线终验

```bash
python3 03_batch_blast_ncbi.py \
    --input candidates.faa \
    --output nr_results \
    --database nr \
    --evalue 1e-5 \
    --hitlist 100
```

- 参数：
  | 参数 | 默认 | 说明 |
  |---|---|---|
  | `--input` | 必填 | 候选蛋白 FASTA |
  | `--output` | 必填 | 输出目录（每条一个 XML） |
  | `--database` | nr | BLAST 数据库 |
  | `--evalue` | 0.001 | E-value 阈值 |
  | `--hitlist` | 500 | 最大返回 hit 数 |
- 逻辑：逐条提交在线 BLASTp → 轮询 RID → 完成后下载 XML
- 输出：`<输出目录>/<蛋白ID>.xml`（XML 含全部 hit 与比对序列）+ `<输出目录>/rids.tsv`（提交记录）
- 特性：
  - **断点续跑**（已完成的蛋白跳过；中断后重跑自动续接 RID）
  - **网络失败自动重试**（每次操作重试 3 次）
  - 支持 NCBI API Key（提升请求频率 3 → 10 次/秒）

---

## 四、判定标准（本流程采用的阈值）

| 环节 | 阈值 | 依据 |
|---|---|---|
| 本地互比显著性 | E-value < 1e-3 | BLAST 常规显著阈值 |
| 分级 A/B/C | identity 30% / 40% | 30% 为远缘同源上限；40% 为种间保守区判定线 |
| 在线终验 | E-value < 1e-5 | nr 库庞大，需更严阈值控制假阳性 |
| 全长同源判定 | identity ≥ 40% 且 query coverage ≥ 50% | 排除片段相似，只计全长同源 |

---

## 五、注意事项

1. **运行顺序**：严格按照 ① → ② → ③ → ④ 执行，后一步依赖前一步输出
2. **断点续跑**：所有脚本支持中断后重跑（自动跳过已完成部分），适合长时任务
3. **磁盘空间**：② 步骤会为每个菌临时建库（约数 GB），脚本结束后自动清理
4. **在线 BLAST 时长**：nr 库单条查询约 1–3 分钟（受 NCBI 队列负载影响）；建议分批提交
5. **长时间运行**：建议在 `tmux`/`screen` 会话中执行，避免终端断开中断任务
   ```bash
   tmux new -s screen_run
   # 执行脚本...
   # Ctrl+B 然后 D 脱离会话；tmux attach -t screen_run 可重新连接
   ```
6. **可复现性**：下载脚本固定了 RefSeq 组装编号（GCF_xxx），保证每次运行获得相同的基因组版本

---

## 六、文件清单

| 文件 | 说明 |
|---|---|
| `01_download_proteomes.sh` | 下载 32 菌 RefSeq 蛋白组 |
| `02_blast_specificity_v2.py` | 32 菌互比（本地 BLASTp） |
| `03_batch_blast_ncbi.py` | NCBI nr 库在线 BLASTp 终验 |
| `04_rank_specificity_graded.py` | 互比结果分级（S/A/B/C） |
| `README.md` | 本文档 |
