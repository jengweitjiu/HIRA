# Claude Code × Breakthrough Strategy 實戰指南

## 核心概念：Claude Code 如何加速你的工作

Claude Code 是終端機裡的 agentic coding assistant——它讀你的整個 codebase、編輯檔案、執行指令、迭代修 bug，全部用自然語言驅動。對你的 breakthrough strategy 而言，它能：

1. **自動寫分析 pipeline**：描述你要做什麼（如「run TOPPLE on CIMA Table S5」），它會生成完整 Python 腳本
2. **跨多檔案協調**：同時管理 TOPPLE、SICAI、DGSA 的程式碼
3. **Debug 迴圈**：跑完 script → 遇到 error → 自動修復 → 再跑，不需你介入
4. **Git 整合**：自動 commit、push 到 GitHub (jengweitjiu/TOPPLE 等)

---

## Step 0: 安裝 Claude Code

```bash
# 需要 Node.js 18+
# macOS
brew install node
# 或用 nvm
nvm install 20

# 安裝 Claude Code
npm install -g @anthropic-ai/claude-code

# 進入你的專案目錄
cd ~/CIMA_Breakthrough
claude
```

需要 Claude Pro、Max、Team 或 Enterprise 訂閱，或者 Anthropic Console API 帳號。

---

## Step 1: 建立專案結構

```bash
mkdir -p ~/CIMA_Breakthrough/{data,scripts,results,figures}

# 把 CIMA supplementary tables 放進 data/
cp science_adt3130_table_s*.xlsx ~/CIMA_Breakthrough/data/
cp science_adt3130_table_s*.csv  ~/CIMA_Breakthrough/data/
```

---

## Step 2: 建立 CLAUDE.md（關鍵步驟）

Claude Code 每次啟動 session 會讀 `CLAUDE.md`，它等於你的「專案記憶」。這是讓 Claude Code 理解你整個研究脈絡的關鍵：

```markdown
# CLAUDE.md — CIMA × Methods Breakthrough Project

## Project Goal
Apply TOPPLE, SICAI, DGSA, IPA frameworks to CIMA (Yin et al., Science 2026)
published supplementary tables. Generate 2-figure proof of concept for
Genome Biology / NatComms submission.

## Key Constraint
- NO raw data access (OMIX007351 requires Chinese institutional access)
- ALL analyses use published supplementary tables ONLY
- Environment: Colab Pro+ (T4 GPU, <64GB RAM) or local Python 3.10+

## Data Files (in data/)
- Table S5: 203 eRegulons × 61 cell types (AUC, RSS, age_corr, sex_diff)
  → 12,383 entries → TOPPLE input
- Table S6: 223,405 lead xQTLs (eQTL + caQTL) → DGSA input
- Table S8: 4,692 eQTL sharing pairs (π₁, r_b) → SICAI input
- Table S15: 2,085 SMR pleiotropic associations → disease correlation

## Methods Reference
- TOPPLE: Regulon stability landscape via perturbation analysis
  Key finding from psoriasis: STAT5B = #1 stabilizer, redistribution index 0.03→0.61
  GitHub: jengweitjiu/TOPPLE (v0.2.2, pySCENIC v2, GSE173706)
- SICAI: Coupling complexity as biomarker
  Key finding: coupling complexity–PASI ρ=0.480, P=0.020, n=23
- DGSA: Geometric decomposition of feature ablation / stability analysis
  GitHub: jengweitjiu/DGSA-stability
- IPA: Perturbation logic for regulatory inference (kill experiment)
  GitHub: jengweitjiu/IPA-kill-experiment

## Coding Standards
- Python 3.10+, scanpy ≥1.9.6, pandas, numpy, matplotlib, seaborn
- Use harmonypy directly (NOT sc.external.pp.harmony — shape bug)
- HVG: seurat_v3 flavor, need scikit-misc
- Figures: publication quality, 300 DPI, PDF + PNG
- Single-author style: all code self-contained, reproducible

## Current Manuscript Status
- SICAI → Genome Biology (submission ready)
- TOPPLE → Nature Methods (manuscript complete, needs refs + submit)
- DGSA → Bioinformatics Advances (ACCEPT TRANSFER — 30-day expiry from 7-Mar)
- STRATA → Communications Medicine (transfer pending)
- STRATA-Atlas → Science Advances (transfer pending)

## Build/Run Commands
pip install scanpy harmonypy scikit-misc pyscenic pandas openpyxl matplotlib seaborn
python scripts/01_topple_cima.py
python scripts/02_sicai_cima.py
python scripts/03_dgsa_cima.py
```

---

## Step 3: 用 Claude Code 執行 72-Hour Plan

### Day 1 — TOPPLE on CIMA Table S5

在 terminal 中啟動 Claude Code 後，直接對話：

```
> Read data/science_adt3130_table_s5.xlsx and extract the 203 regulons × 61 
  cell types AUC matrix. Then implement TOPPLE stability analysis:
  1. Build regulon activity matrix from AUC scores
  2. Compute perturbation response for each regulon removal
  3. Calculate redistribution index per regulon per cell type
  4. Identify universal stabilizers and destabilizers across immune lineages
  5. Check if STAT5B emerges as top stabilizer (matches psoriasis finding)
  6. Generate publication figure: stability landscape heatmap
  Save script to scripts/01_topple_cima.py and figure to figures/
```

Claude Code 會：
- 讀取 Excel 檔案
- 寫完整 Python script
- 執行它
- 遇到 error 自動修復
- 產出 figure

### Day 2 — SICAI on CIMA Table S8

```
> Load data/science_adt3130_table_s8 (4,692 eQTL sharing pairs with r_b values).
  Build the 69×69 cell type r_b correlation matrix.
  Compute SICAI coupling complexity per cell type.
  Then load Table S15 (2,085 SMR pleiotropic associations) and count 
  disease associations per cell type.
  Test correlation: does coupling complexity predict disease pleiotropy?
  Generate scatter plot with Spearman rho and P-value.
  Save to scripts/02_sicai_cima.py and figures/
```

### Day 3 — Combine & Draft

```
> Combine the TOPPLE and SICAI results into a 2-figure proof of concept:
  Fig 1: TOPPLE stability landscape across CIMA's 61 cell types
  Fig 2: Coupling complexity vs. disease association count
  Create a multi-panel figure (2 rows) suitable for Genome Biology.
  Also draft a 500-word framing paragraph explaining why these architectural
  properties were invisible in CIMA's original analysis.
  Save figure to figures/fig_combined_proof.pdf
  Save text to results/genome_biology_framing.md
```

---

## Step 4: 進階 Claude Code 技巧

### 用 Custom Commands 包裝重複工作流程

在 `.claude/commands/` 建立自定義指令：

```markdown
# .claude/commands/run-topple.md
Run TOPPLE stability analysis on the specified input file.
Load the regulon × cell type AUC matrix, compute perturbation response,
calculate redistribution index, identify stabilizers/destabilizers,
and generate publication-quality heatmap.
```

之後只要輸入 `/run-topple data/new_dataset.xlsx`

### 用 Sub-agents 平行處理

```
> Use subagents to parallelize:
  Agent 1: Run TOPPLE on Table S5 (stability analysis)
  Agent 2: Run SICAI on Table S8 (coupling complexity)  
  Agent 3: Run DGSA on Table S6 (geometric decomposition)
  Each agent saves results to results/ directory.
  Then merge all three into a unified summary.
```

### 搭配 MCP Servers 整合外部工具

```bash
# 在 ~/.claude/settings.json 加入 MCP servers
{
  "mcpServers": {
    "github": {
      "command": "npx",
      "args": ["-y", "@modelcontextprotocol/server-github"],
      "env": { "GITHUB_TOKEN": "your_token" }
    }
  }
}
```

這樣 Claude Code 可以直接：
- 讀取 jengweitjiu/TOPPLE repo 的 issues
- 提交 PR
- Push code 到 GitHub

---

## Step 5: Colab 整合策略

因為你主要用 Colab Pro+，最佳工作流是：

### 方案 A：Claude Code 寫 → Colab 跑（推薦）
```
Local (Claude Code)          Colab Pro+
┌──────────────────┐        ┌──────────────────┐
│ 寫 scripts       │───────→│ 上傳 & 執行      │
│ Debug logic      │        │ GPU/大記憶體運算  │
│ Git push         │        │ 結果下載         │
│ 整理 figures     │←───────│                  │
└──────────────────┘        └──────────────────┘
```

1. Claude Code 在本機寫好完整 .py scripts
2. Push 到 GitHub
3. Colab `!git clone` 或 `!pip install` 後直接跑
4. 結果下載回本機，Claude Code 幫你整理 figures

### 方案 B：Claude Code 直接跑（如果 table 分析不需 GPU）
Table S5/S8/S15 的分析都是 pandas/numpy 運算，不需 GPU：
- Claude Code 直接在本機跑完
- 只有需要 Harmony/pySCENIC 的大規模運算才上 Colab

---

## Step 6: 具體 Prompt 範本庫

### Prompt 1: 驗證 TOPPLE 跨 atlas 一致性
```
Compare TOPPLE results from CIMA (Table S5, 61 cell types, 203 regulons) 
with our psoriasis results (GSE173706, 83,352 cells, 1,492 regulons).
Focus on: does STAT5B rank as top stabilizer in both?
Generate a Venn diagram of top-20 stabilizers in each dataset.
```

### Prompt 2: DGSA 幾何分解
```
Load Table S6 (223,405 xQTLs). For each gene with xQTLs in multiple 
cell types, decompose the effect size vector geometrically.
Compute the non-additivity index across cell types.
Identify genes where the geometric structure deviates most from 
simple additive models. These are candidates for regulatory architecture.
```

### Prompt 3: 一鍵產出 submission package
```
Package the SICAI Genome Biology submission:
1. Check Fig3c-d labels are correct
2. Verify Fig5 panel map
3. Check Reference 11
4. Estimate APC (~$3,290)
5. Create SICAI_GenomeBiology_Submission.zip with all figures, 
   supplementary tables, and manuscript
```

### Prompt 4: HIRA 統一框架
```
Read the existing manuscripts for TOPPLE, DGSA, SICAI, IPA, STRATA, 
and STRATA-Atlas. Create a unified HIRA framework document that shows
how the five layers connect:
TOPPLE → DGSA → SICAI → STRATA → IPA+STRATA-Atlas
Map each layer to specific CIMA tables and Thompson et al. data.
Output as a Mermaid diagram + narrative.
```

---

## 重要提醒

| 優先順序 | 行動 | 工具 |
|---------|------|------|
| **🔴 立即** | Accept DGSA → Bioinformatics Advances transfer | 手動（期限 30 天，從 3/7 起算） |
| **Day 1** | TOPPLE on CIMA Table S5 | Claude Code |
| **Day 2** | SICAI on CIMA Table S8 + S15 | Claude Code |
| **Day 3** | Combined figure + Genome Biology framing | Claude Code |
| **本週** | SICAI pre-submit checklist (Fig3c-d, Fig5, Ref11) | Claude Code |
| **本週** | TOPPLE: add references, Zenodo deposit, GitHub push, submit | Claude Code + 手動 |

---

## Claude Code vs Claude.ai (這裡) 的分工

| 任務類型 | 用 Claude Code | 用 Claude.ai |
|---------|---------------|-------------|
| 寫 Python/R scripts | ✅ 自動執行+debug | ❌ 只能顯示程式碼 |
| 跨多檔案編輯 | ✅ 一次改多個 .py | ❌ 一次一個檔案 |
| Git 操作 | ✅ commit/push/PR | ❌ 無法 |
| 讀大型 Excel tables | ✅ 直接 pandas 操作 | ⚠️ 限制較大 |
| 策略規劃 & 寫作 | ⚠️ 可以但非強項 | ✅ 長文脈絡更好 |
| Manuscript 寫作/潤稿 | ⚠️ 適合結構化任務 | ✅ 更適合 |
| Figure 美化 | ✅ matplotlib 迭代 | ⚠️ 只能給建議 |
| Colab notebook 生成 | ✅ 直接寫 .ipynb | ⚠️ 有限 |
