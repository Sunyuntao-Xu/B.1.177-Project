# 🧬 B.1.177-Project: SARS-CoV-2 Alignment & Phylodynamic Analysis Pipeline

This repository contains a complete alignment preparation and modeling pipeline for the **B.1.177 lineage of SARS-CoV-2**, including:

- Reference construction and sequence cleaning
- Key mutation inspection and site masking
- Finalized alignment for phylodynamic modeling
- Case smoothing and migration modeling using **PhydynR**

---

## 📁 Folder Structure

### `alignment/`
Scripts for **sequence alignment**, **reference preparation**, **mutation inspection**, and **site masking**.

| Script                                 | Description |
|----------------------------------------|-------------|
| `B.1.177_alignment_cleanup_pipeline.sh` | Main pipeline script for processing: merges references, extracts B.1.177 sequences, trims alignments, and performs initial cleanup. |
| `remove_gap_only_columns.py`           | Python script to remove columns in the alignment that consist only of gaps across all sequences. |
| `B.1.177_key_mutation_retrive.sh`      | Extracts nucleotide states at selected key mutation positions for inspection and tracking. |
| `B.1.177_sites_check_and_masking.sh`   | Summarizes base frequencies at predefined mutation sites and masks specific sites using `seqkit mutate`. |
| `filter_B.1.177_by_mutations.sh`       | Filters the alignment to include only sequences carrying defined mutations. |
| `split_alignment.sh`                   | Optional utility to split large FASTA files into smaller chunks for easier handling or parallel processing. |
| `subsampling_and_modelfinding.sh`      | Performs alignment subsampling and selects the best substitution model using IQ-TREE. |
| `B.1.177_tree_building.sh`             | Constructs maximum likelihood trees using IQ-TREE, with support for bootstrapping and model selection. |


---

## 📜 Main Analysis Scripts

### `PhydynR_Test_Code_and_UShER_Selection.R`
- Validates UShER reference trees
- Aligns and dates reference sequences
- Sets up initial `phydynR` test runs

### `Project_phydynR_for_migration_modelling.R`
- Main phylodynamic modeling script
- Implements **SIR model with migration** across geographic demes
- Incorporates **time-varying transmission rate** `β(t)`
- Incorporates **time-varying migration rate** `m(t)`
- Calibrated using alignment and smoothed case data

### `GISAID.R`
- Parses full GISAID metadata
- Selects B.1.177 sequences by ID list
- Optionally filters by country, date, or lineage

### `B.1.177_case_smoothing.R`
- Applies **Generalized Additive Models (GAMs)** to UKHSA/COG-UK data
- Estimates B.1.177 growth rate and inflection points
- Provides time windows for time-stratified modeling in `phydynR`

---
## 📊 Required Input Files

| File | Description |
|------|-------------|
| `msaCodon_0201_fixed.fasta` | Full codon-aware alignment from GISAID |
| `GISAID_selected_sequences.txt` | List of selected B.1.177 sequence headers |
| `Wuhan-Hu-1.fasta` + 11 UShER reference FASTAs | Reference sequences for anchoring alignment |
| UKHSA/COG-UK CSV case data | For smoothing and time-series modeling |

---

## ⚙️ Tools & Dependencies

- [`MAFFT`](https://mafft.cbrc.jp/alignment/software/): reference and sample alignment
- [`seqkit`](https://bioinf.shenwei.me/seqkit/): FASTA filtering, mutation masking
- [`Biopython`](https://biopython.org/): gap-only column removal
- [`R` & `phydynR`](https://github.com/emvolz-phylodynamics/phydynR): migration + phylodynamic modeling
- [`PowerShell`](https://learn.microsoft.com/en-us/powershell/): Windows-based sequence inspection and mutation grouping
- [`trimAl`](https://vicfero.github.io/trimal/):automated removal of highly gappy or poorly conserved alignment columns

---

## 🔁 Pipeline Summary

1. **Reference Preparation**
    - Merge Wuhan-Hu-1 and UShER sequences
    - Align with MAFFT
    - Manually trim ends if needed

2. **GISAID Sequence Extraction**
    - Select B.1.177 sequences via `seqkit grep`
    - Align to references using `mafft --add --keeplength`

3. **Alignment Cleaning**
    - Remove gap-only columns (`remove_gap_only_columns.py`)
    - Mask problematic or hypervariable sites (`seqkit mutate`)

4. **Mutation Inspection**
    - Extract base state at key SNPs
    - Group and summarize by base using PowerShell

5. **Phylodynamic Modeling**
    - Smooth case data using GAM
    - Construct SIR + migration model using `phydynR`
    - Time-stratified modeling of B.1.177 spread

---

## 📌 Example Usage

```bash
bash alignment/B.1.177_alignment_cleanup_pipeline.sh
python alignment/remove_gap_only_columns.py
bash alignment/B.1.177_sites_check_and_masking.sh
Rscript phylodynamics/Project_phydynR_for_migration_modelling.R
