# Merge Wuhan-Hu-1 and 11 UShER references into one FASTA file
cat \
  "E:/Reference Alignments/Wuhan-Hu-1.fasta" \
  "E:/Reference Alignments/MT991088.1.fasta" \
  "E:/Reference Alignments/MW056032.1.fasta" \
  "E:/Reference Alignments/MW404673.1.fasta" \
  "E:/Reference Alignments/OA970147.1.fasta" \
  "E:/Reference Alignments/ON194679.1.fasta" \
  "E:/Reference Alignments/ON282567.1.fasta" \
  "E:/Reference Alignments/ON302325.1.fasta" \
  "E:/Reference Alignments/ON302991.1.fasta" \
  "E:/Reference Alignments/OU076971.1.fasta" \
  "E:/Reference Alignments/LR881876.2.fasta" \
  "E:/Reference Alignments/LR882266.2.fasta" \
  > "E:/Reference Alignments/merged_reference_sequences.fasta"

# Align merged reference sequences using MAFFT
mafft --auto "E:/Reference Alignments/merged_reference_sequences.fasta" > "E:/Reference Alignments/merged_reference_sequences_aligned.fasta"


# Extract B.1.177 sequences using SeqKit
seqkit grep -f "E:/GISAID_selected_sequences.txt" "E:/msaCodon_0201_fixed.fasta" > "E:/msaCodon_0201_B.1.177.fasta"

# Use trimAl to remove highly gappy columns (>10% gaps) from B.1.177 sequences while preserving at least 95% of the alignment,
# ensuring cleaner alignment and more consistent structure before aligning to the reference
cd C:\trimAl
.\trimal.exe -in "E:\msaCodon_0201_B.1.177.cleaned.fasta" -out "E:\msaCodon_0201_B.1.177.final_trimmed.fasta" -gt 0.9 -cons 95

# Use seqkit to clean out all gaps
seqkit replace -p " " -r "_" \
  -o E:/msaCodon_0201_B.1.177.final_trimmed.fasta \
  E:\msaCodon_0201_B.1.177.final_trimmed_nospace.fasta


# Align B.1.177 sequences to masked reference using MAFFT --add
mafft --keeplength --add "E:/msaCodon_0201_B.1.177.cleaned.fasta" \
      "E:/Reference Alignments/merged_reference_sequences_aligned_masked.fasta" \
      > "E:/aligned_B.1.177_to_reference.fasta"

