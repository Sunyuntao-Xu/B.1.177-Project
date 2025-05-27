# Position 22227 → Expect T
$bad1 = Import-Csv "E:/pos_22227_check.tsv" -Delimiter "`t" -Header "ID", "Base" |
    Where-Object { $_.Base -ne "T" } | Select-Object -ExpandProperty ID

# Position 23403 → Expect G
$bad2 = Import-Csv "E:/pos_23403_check.tsv" -Delimiter "`t" -Header "ID", "Base" |
    Where-Object { $_.Base -ne "G" } | Select-Object -ExpandProperty ID

# Position 28932 → Expect T
$bad3 = Import-Csv "E:/pos_28932_check.tsv" -Delimiter "`t" -Header "ID", "Base" |
    Where-Object { $_.Base -ne "T" } | Select-Object -ExpandProperty ID

# Position 29645 → Expect T
$bad4 = Import-Csv "E:/pos_29645_check.tsv" -Delimiter "`t" -Header "ID", "Base" |
    Where-Object { $_.Base -ne "T" } | Select-Object -ExpandProperty ID

# Combine and write cleanly
$allBad = $bad1 + $bad2 + $bad3 + $bad4 | Sort-Object -Unique
$allBad | Set-Content -Encoding UTF8 "E:/B.1.177_exclude_ids.txt"



# Fix the GISAID ID list by replacing spaces with underscores
Get-Content "E:\GISAID_selected_sequences.txt" |
    ForEach-Object { $_.Trim() -replace " ", "_" } |
    Where-Object { $_ -ne "" } |
    Set-Content "E:\GISAID_selected_sequences_nospace.txt"


Get-Content "E:\B.1.177_exclude_ids.txt" |
    ForEach-Object { $_.Trim() -replace " ", "_" } |
    Where-Object { $_ -ne "" } |
    Set-Content "E:\B.1.177_exclude_ids_nospace.txt"

seqkit grep -f "E:\B.1.177_exclude_ids_nospace.txt" -v "E:\masked_output.fasta" -o "E:\key_site_seqkit_filtered_B.1.177.fasta"

# REmove duplicate IDs
seqkit grep -f E:/iqtree_duplicate_ids.txt -v E:/key_site_seqkit_filtered_B.1.177.fasta -o E:/key_site_seqkit_filtered_B.1.177.dedup.fasta

# Identify identical sequences
seqkit rmdup -s -i -o NUL -d duplicated.fa.gz -D E:/identical.detail.txt E:/key_site_seqkit_filtered_B.1.177.dedup.fasta






