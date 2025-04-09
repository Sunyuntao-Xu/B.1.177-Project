# B.1.177-Project

This repository contains the full analysis pipeline and related scripts for the **B.1.177 SARS-CoV-2 lineage** project, focusing on:

- Reference alignment and GISAID sequence preparation
- Case smoothing and incidence estimation
- PhydynR-based migration and phylodynamic modeling
- UShER reference validation
- Data cleaning and preparation

---

## 📁 Folder Overview

### `alignment/`
Contains all code related to the alignment workflow, including:

- Merging Wuhan-Hu-1 with 11 UShER reference sequences
- Aligning references using MAFFT
- Manual masking of alignment ends
- Selecting B.1.177 sequences from GISAID using SeqKit
- Aligning B.1.177 sequences to masked references using `mafft --add --keeplength`
- Python script to remove gap-only columns from final alignment

**Key scripts:**
- `alignment-workflow.sh`: Full pipeline summary in shell script form
- `Delete gap only columns for alignments.py`: Removes columns that are gaps in all sequences

---

## 📜 Script Overview

### `PhydynR Test Code and UShER Selection.R`
- Tests initial setup of PhydynR
- Performs UShER reference validation and tree dating for reference selection

### `Project phydynR for migration modelling.R`
- Main PhydynR modeling script
- Implements **SIR compartment model with migration**
- Includes time-stratified beta(t) function and deme migration matrices

### `GISAID.R`
- Parses and filters GISAID metadata
- Matches sequence IDs to external lists (e.g. for B.1.177 selection)

### `B.1.177 case smoothing.R`
- Smooths case incidence data from COG-UK and UKHSA
- Uses GAM (Generalized Additive Models) to estimate growth rates
- Identifies key inflection points in B.1.177 spread for time-stratified modeling

---

## 📊 Data Requirements

**Required Input Files:**
- GISAID full alignment: `msaCodon_0201_fixed.fasta`
- GISAID metadata selection list: `GISAID_selected_sequences.txt`
- Reference sequences (FASTA format): Wuhan-Hu-1 + 11 UShER isolates
- UKHSA/COG-UK case data (CSV format) for smoothing

---

## ⚙️ Tools Used

- [MAFFT](https://mafft.cbrc.jp/alignment/software/): multiple sequence alignment
- [SeqKit](https://bioinf.shenwei.me/seqkit/): FASTA/FASTQ filtering
- [Biopython](https://biopython.org/): gap-removal script
- [R](https://www.r-project.org/): smoothing, phylodynamic modeling
- [PhydynR](https://cran.r-project.org/web/packages/phydynR/index.html): compartmental modeling from phylogenies
