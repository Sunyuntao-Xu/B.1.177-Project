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
seqkit grep --by-name -f E:/GISAID_selected_sequences.txt E:/msaCodon_0201_fixed.fasta > E:\B.1.177_extracted.fasta

# Alternative B.1.177 sequences using bash command

# Set file paths
$fastaFile = "E:\msaCodon_0201_fixed.fasta"
$nameListFile = "E:\GISAID_selected_sequences.txt"
$outputFile = "E:\B.1.177_extracted.fasta"

# Read target headers into a hashset for fast lookup
$targetNames = Get-Content $nameListFile | ForEach-Object { $_.Trim() } | Where-Object { $_ -ne "" }
$targetSet = @{}
foreach ($name in $targetNames) {
    $targetSet[$name] = $true
}

# Read FASTA and extract matching sequences
$writeBlock = $false
Get-Content $fastaFile | ForEach-Object {
    if ($_ -like ">*") {
        $header = $_.Substring(1).Trim()
        if ($targetSet.ContainsKey($header)) {
            $writeBlock = $true
            Add-Content -Path $outputFile -Value $_
        }
        else {
            $writeBlock = $false
        }
    }
    elseif ($writeBlock) {
        Add-Content -Path $outputFile -Value $_
    }
}

# Use trimAl to remove highly gappy columns (>10% gaps) from B.1.177 sequences while preserving at least 95% of the alignment,
# ensuring cleaner alignment and more consistent structure before aligning to the reference
cd C:\trimAl
.\trimal.exe -in "E:\msaCodon_0201_B.1.177.cleaned.fasta" -out "E:\msaCodon_0201_B.1.177.final_trimmed.fasta" -gt 0.9 -cons 95



# Align B.1.177 sequences to masked reference using MAFFT --add
mafft --keeplength --add "E:/msaCodon_0201_B.1.177.cleaned.fasta" \
      "E:/Reference Alignments/merged_reference_sequences_aligned_masked.fasta" \
      > "E:/aligned_B.1.177_to_reference.fasta"

