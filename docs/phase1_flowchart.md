# Phase 1 flowchart

How Phase 1 (two-sample MR nomination) runs from inputs to the final tiered
shortlist, with the decision points at each step. Script numbers refer to
`scripts/R/NN_*.R`; column definitions are in
[`phase1_output_spec.md`](phase1_output_spec.md). Counts are for the primary
analysis (controls mQTL meta × Nalls 2019 GCST009325).

## 1. Overview

```mermaid
flowchart TD
    subgraph inputs [Inputs]
        MQTL["cis-mQTL summary stats<br/>PD-QTLs: PEG1 + PPMI<br/>controls meta (primary), cases meta,<br/>single-cohort controls"]
        GWAS["PD GWAS, Nalls 2019<br/>GCST009325 (primary)<br/>GCST009324 no-proxy (sensitivity)"]
        REF["1000G phase 3 EUR<br/>LD reference"]
        IND["PEG / PPMI individual-level<br/>genotypes, methylation, covariates<br/>(local only)"]
    end

    S1["Step 1: exposure candidates<br/>cis pairs p < 5e-8, probe-SNP filter"]
    S2["Step 2: harmonise to GWAS<br/>alleles, strand, palindromes"]
    S3["Step 3: clump per CpG<br/>r² < 0.001"]
    S4["Step 4: MR per CpG<br/>Wald / IVW, Steiger, BH FDR"]
    HIT{"FDR < 0.05?"}
    STOP1(["not nominated<br/>12,734 CpGs"])
    S5["Step 5: regional mQTL refit<br/>±1 Mb, untruncated"]
    S6["Step 6: coloc.abf<br/>one causal variant"]
    S7["Step 7: coloc-SuSiE<br/>several causal variants"]
    S8["Step 8: in-house GWAS<br/>PPMI, PEG"]
    X["Crossed replication<br/>PEG mQTL × PPMI GWAS<br/>PPMI mQTL × PEG GWAS"]
    S10["Step 10: MAPT H1/H2<br/>conditional analyses"]
    S11["Step 11: shared instruments<br/>group CpGs by instrument LD"]
    S12["Step 12: probe-body SNPs<br/>flag instruments in LD<br/>with a SNP in the probe"]
    S9["Step 9: results table<br/>evidence tier A / B / C"]
    OUT(["Nominations<br/>44 CpGs in 22 instrument groups<br/>tier A: 5 groups, B: 5, C: 12"])

    MQTL --> S1 --> S2
    GWAS --> S2 --> S3
    REF --> S3 --> S4 --> HIT
    HIT -- "no" --> STOP1
    HIT -- "yes: 44 CpGs" --> S5
    IND --> S5 --> S6 --> S7 --> S9
    IND --> S8 --> X --> S9
    S5 --> S10 --> S9
    REF --> S10
    S5 --> S11 --> S9
    S5 --> S12 --> S9
    REF --> S12
    REF --> S11
    S9 --> OUT
```

Mandatory run order: 1 → 2 → 3 → 4 → 5 → 6 → 7 → 8 (steps 1–4 again for the
crossed pairings) → 10 → 11 → 12 → 9.

## 2. Building instruments (steps 1–3)

Applied to every (SNP, CpG) pair in each exposure dataset.

```mermaid
flowchart TD
    P["cis-mQTL pair<br/>SNP within 1 Mb of CpG"]
    D1{"mQTL p < 5e-8?"}
    D2{"SNP more than 50 bp<br/>from the CpG?"}
    D3{"Meta exposure:<br/>same REF/ALT in<br/>both cohorts?"}
    D4{"SNP at the same<br/>chr:pos in the GWAS?"}
    D5{"Alleles match<br/>(direct, swapped,<br/>or strand-flipped)?"}
    D6{"Palindromic<br/>(A/T or C/G)?"}
    D7{"Exposure MAF < 0.42<br/>(frequency resolves<br/>the strand)?"}
    D8{"abs(eaf_x − eaf_y) ≤ 0.2?"}
    H["harmonised SNP<br/>(GWAS effect re-signed<br/>to the ALT allele)"]
    C1{"Lead SNP for this CpG<br/>in 1000G EUR?"}
    C2["Clump within the CpG:<br/>keep lead, drop SNPs with<br/>r² ≥ 0.001 to a kept SNP, repeat"]
    C3["use the strongest SNP<br/>that is in the reference<br/>(ld_ref = lead_not_in_ref)"]
    C4["no candidate in reference:<br/>keep top SNP unclumped<br/>(ld_ref = missing)"]
    IV(["instruments for the CpG<br/>12,778 CpGs; 99.6% have exactly 1"])
    X(["dropped"])

    P --> D1
    D1 -- no --> X
    D1 -- yes --> D2
    D2 -- "no: probe artefact" --> X
    D2 -- yes --> D3
    D3 -- "no: multiallelic mismatch" --> X
    D3 -- "yes / single cohort" --> D4
    D4 -- "no: missing from GWAS" --> X
    D4 -- yes --> D5
    D5 -- no --> X
    D5 -- yes --> D6
    D6 -- no --> D8
    D6 -- yes --> D7
    D7 -- "no: ambiguous" --> X
    D7 -- yes --> D8
    D8 -- "no: frequency discordant" --> X
    D8 -- yes --> H --> C1
    C1 -- yes --> C2 --> IV
    C1 -- "no, others are" --> C3 --> C2
    C1 -- "none in reference" --> C4 --> IV
```

About 2% of exposure SNPs are lost at the GWAS steps (`harmonised/summary.tsv`).
The GWAS p-value is never used to choose instruments.

## 3. MR estimate per CpG (step 4)

```mermaid
flowchart TD
    I["instruments for one CpG"]
    N{"How many<br/>instruments?"}
    W["Wald ratio<br/>b = beta_y / beta_x"]
    IVW["IVW, multiplicative<br/>random effects + Cochran's Q"]
    R["+ weighted median<br/>+ MR-Egger (intercept test)"]
    ST["Steiger direction check<br/>r²(methylation) > r²(PD)?"]
    F["Benjamini–Hochberg FDR<br/>across all CpGs in this<br/>exposure × outcome"]
    Q{"FDR < 0.05?"}
    H(["MR hit → steps 5–11"])
    NO(["not nominated"])

    I --> N
    N -- "1 (12,730 CpGs)" --> W --> ST
    N -- "2 (48 CpGs)" --> IVW --> ST
    N -- "≥ 3 (0 in primary)" --> IVW
    IVW -. "≥ 3 only" .-> R --> ST
    ST --> F --> Q
    Q -- yes --> H
    Q -- no --> NO
```

The same steps 1–4 run for every exposure × outcome pairing: the sensitivity
exposures (cases meta, PEG1 controls, PPMI controls), the no-proxy GWAS, and the
crossed in-house pairings.

## 4. Evidence tier (steps 6, 7, 10 → 9)

MR with one instrument can't tell "the same variant drives methylation and PD"
from "two different variants in LD". Colocalization decides that, and sets the tier.

```mermaid
flowchart TD
    H["MR hit (44 CpGs)"]
    S{"coloc-SuSiE ran and<br/>found a credible set<br/>for both traits?"}
    SH4{"SuSiE best H4"}
    A{"coloc.abf<br/>(fallback)"}
    TA["tier A"]
    TB["tier B"]
    TC["tier C"]
    M{"In the MAPT<br/>17q21.31 region?"}
    MC{"H2-conditional<br/>coloc H4 ≥ 0.8?"}
    FA(["final tier A"])
    FB(["final tier B"])
    FC(["final tier C"])

    H --> S
    S -- yes --> SH4
    SH4 -- "≥ 0.8" --> TA
    SH4 -- "0.5–0.8" --> TB
    SH4 -- "otherwise (H3 favoured)" --> TC
    S -- "no credible set /<br/>not run" --> A
    A -- "H4 ≥ 0.8" --> TA
    A -- "H4 0.5–0.8" --> TB
    A -- "otherwise" --> TC
    TA --> M
    TB --> M
    TC --> M
    M -- no --> FA
    M -- no --> FB
    M -- no --> FC
    M -- yes --> MC
    MC -- "yes (0 CpGs)" --> FA
    MC -- "no (12 CpGs)" --> FC
```

Not shown because they don't change the tier: coloc-SuSiE with the no-proxy GWAS,
crossed replication, the sensitivity exposures, and the probe-body SNP flag
(step 12, `probe_snp_ld`). They are reported as columns
next to the tier.

## 5. MAPT H1/H2 classification (step 10)

The 17q21.31 inversion makes the H1/H2 haplotype a single ~1.5 Mb LD block. The H2
tag (rs8070723, G allele) is itself a genome-wide mQTL for 45 CpGs.

```mermaid
flowchart TD
    C["MAPT-region MR hit"]
    R{"Instrument in 1000G EUR?"}
    R2{"r² with H2 tag ≥ 0.5?"}
    CM{"Any mQTL p < 5e-8 left<br/>after adding H2 dosage<br/>to the mQTL model?"}
    H2H(["H2_haplotype<br/>6 CpGs: the MR estimate is<br/>the haplotype effect, rescaled"])
    H2R["H2_tagging_residual<br/>2 CpGs"]
    NH["not_H2<br/>4 CpGs"]
    LN(["lead_not_in_ref"])
    CG["H2-conditional GWAS<br/>(approx., 1000G EUR LD)<br/>+ conditional MR + coloc"]
    D{"conditional<br/>coloc H4 ≥ 0.8?"}
    KEEP(["CpG-specific signal:<br/>keep coloc tier"])
    TC(["tier C<br/>(all 12 primary MAPT hits)"])

    C --> R
    R -- no --> LN
    R -- yes --> R2
    R2 -- yes --> CM
    CM -- no --> H2H --> TC
    CM -- yes --> H2R --> CG
    R2 -- no --> NH --> CG
    CG --> D
    D -- yes --> KEEP
    D -- "no / no conditional mQTL" --> TC
```

## 6. Shared instruments (step 11)

```mermaid
flowchart TD
    H["44 primary hits,<br/>lead instrument each"]
    L{"Two CpGs' instruments<br/>identical or EUR r² ≥ 0.8?"}
    G["link them<br/>(single linkage: chains merge)"]
    S(["singleton group<br/>16 groups"])
    MG["multi-CpG group<br/>6 groups, 28 CpGs"]
    MV{"independent instruments<br/>(r² < 0.01) ≥ CpGs<br/>in group?"}
    MVMR(["multivariable MR possible<br/>(0 groups)"])
    DESC["describe the group instead:<br/>co-methylation (residual r),<br/>mQTL–mQTL coloc, genes"]
    U(["group = one nomination;<br/>CpGs within it are not<br/>distinguishable by MR"])

    H --> L
    L -- yes --> G --> MG
    L -- "no partner" --> S
    MG --> MV
    MV -- yes --> MVMR
    MV -- "no (all 6)" --> DESC --> U
```
