# CIMA × Your Methods: Breakthrough Strategy
## Date: 2026-03-10

---

## The Core Insight

CIMA (Science 2026) built the map — 10M cells, 73 cell types, 428 individuals.
Your papers (TOPPLE, SICAI, DGSA, IPA) build the compass — mathematical frameworks
that reveal architectural properties invisible to cataloging approaches.

**CIMA's 18 published supplementary tables already contain everything needed
to validate your methods at Science-paper scale. No raw data access required.**

---

## Manuscript Portfolio Status (as of 10-Mar-2026)

| Paper | Target | Status | Next Action |
|-------|--------|--------|-------------|
| SICAI | Genome Biology | Submission ready | Pre-submit checks |
| TOPPLE | Nature Methods | Manuscript complete | Add references, submit |
| STRATA | Commun Med (transfer) | Decision pending | Wait |
| STRATA-Atlas | Science Advances (transfer) | Decision pending | Wait |
| IPA | TBD | Desk rejected NatComms | Consider merge |
| DGSA | Bioinf Advances (transfer) | **Accept transfer (30-day expiry)** | Act now |
| JAAD Reviews ×2 | JAAD Reviews | Submitted | Wait |

**Total desk rejections: 6 across 6 manuscripts. Pattern = framing problem, not science problem.**

---

## The Problem

Methods demonstrated on small/simulated datasets → editors can't judge impact → desk reject.

CIMA with zero methodological novelty → Science because:
- 428 individuals (scale)
- 60 co-authors (institutional credibility)
- Resource framing (community utility)
- Population-first (Chinese immune atlas)

---

## Three Breakthrough Strategies

### Strategy 1: "The Analytical Layer" Paper (HIGHEST PRIORITY)

**Title concept:** "Regulatory stability architecture of the human immune system
revealed by multi-scale geometric analysis of the CIMA atlas"

**Uses CIMA published tables (already downloaded):**

| Your Method | CIMA Table | Analysis |
|-------------|-----------|----------|
| TOPPLE | S5: 203 regulons × 61 cell types (AUC, RSS) | Stability landscape → stabilizers/destabilizers |
| SICAI | S8: 4,692 r_b pairs across 69 cell types | Coupling complexity → disease pleiotropy correlation |
| DGSA | S6: 223,405 xQTLs with effect sizes | Geometric non-additivity across cell types |
| IPA | S5: sex/age perturbation data (21,290 entries) | Natural perturbation → cascade analysis |

**Why it works:**
- Borrows CIMA's scale (428 individuals, 73 cell types) without needing raw data
- "Analysis of Science-published atlas" framing carries credibility
- Shows what CIMA's catalog *means* — interpretive, not duplicative
- Fills a genuine gap: CIMA never analyzes its own data's architecture

**Target:** Genome Biology, Nature Communications, Cell Systems

**Timeline:** Results achievable in 72 hours on Colab Pro+

### Strategy 2: Merge TOPPLE + DGSA + IPA

**Title concept:** "A Geometric Framework for Multi-Scale Regulatory Analysis"

- One comprehensive paper is harder to desk-reject than three thin ones
- TOPPLE (stability) + DGSA (geometry) + IPA (perturbation) = unified framework
- Each component validates the others
- Reviewers see the full picture, not fragments

**Target:** Nature Methods, Genome Biology

### Strategy 3: The Dermatology Anchor

Add to Strategy 1:
1. Methods validated on CIMA immune atlas (scale credibility)
2. Applied to psoriasis spatial transcriptomics (GSE202011, your existing data)
3. Coupling complexity predicts PASI (ρ=0.480, P=0.020, n=23)

Clinician-scientist doing computational biology = extremely rare positioning.
General method + disease application = hard for editors to dismiss.

---

## Immediate 72-Hour Action Plan

### Day 1 (Tuesday):
```
Load CIMA Table S5 in Colab
→ Extract 203 regulons × 61 cell types AUC matrix
→ Run TOPPLE stability analysis
→ Identify universal stabilizers/destabilizers across immune lineages
→ Compare with CIMA's own findings (do they match STAT5B from psoriasis?)
```

### Day 2 (Wednesday):
```
Load CIMA Table S8 in Colab
→ Build 69×69 r_b matrix
→ Compute SICAI coupling complexity per cell type
→ Correlate with number of SMR disease associations from S15
→ Test: does coupling complexity predict disease pleiotropy?
```

### Day 3 (Thursday):
```
Combine results into a 2-figure proof of concept
→ Fig 1: TOPPLE stability landscape across CIMA's 61 cell types
→ Fig 2: Coupling complexity vs. disease association count
→ Draft a 500-word framing for Genome Biology
```

### Also on Day 1:
**Accept the DGSA → Bioinformatics Advances transfer.** Don't let it expire.

---

## The Asymmetry That Works In Your Favor

CIMA invested millions in data generation but used off-the-shelf analysis.
You invest zero in data generation but bring genuinely new analytical concepts.

When you apply your methods to their data, you get:
- Their scale (428 individuals) without their cost
- Your novelty (stability/geometry/coupling) without their 60 co-authors
- A paper that couldn't exist without both contributions
- A framing that editors recognize ("analysis of major atlas")

**CIMA built the telescope. You're the one who looks through it and sees structure.**

---

## Key CIMA Table Reference (for quick access)

| Table | Contents | Size | Your Use |
|-------|----------|------|----------|
| S1 | Metadata, 73 cell type markers | 428 individuals | Cell type definitions |
| S5 | 203 eRegulons: AUC, RSS, age corr, sex diff | 12,383 entries | **TOPPLE input** |
| S6 | 223,405 lead xQTLs (eQTL + caQTL) | 69+42 cell types | **DGSA input** |
| S7 | 13,826 caQTL-eQTL SMR pairs | 40 cell types | Peak→gene regulatory logic |
| S8 | 4,692 eQTL sharing pairs (π₁, r_b) | 69×68 pairs | **SICAI input** |
| S11 | 84 trans-eGenes | 45 cell types | Cross-regulation |
| S12 | 3,126 dynamic eQTLs (B + Monocyte) | LRT P-values | Temporal stability |
| S15 | 2,085 SMR pleiotropic associations | 68 cell types × 68 traits | **Disease correlation** |

All files: already uploaded to Colab as science_adt3130_table_s*.xlsx/csv/zip
