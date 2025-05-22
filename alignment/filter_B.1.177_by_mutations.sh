# Define input and output paths
$input = "E:/masked_output.fasta"
$output = "E:/key_site_filtered_B.1.177.fasta"

# Position 1: 22173 (expect T)
seqkit subseq -r 22173:22173 $input | seqkit fx2tab | where { $_ -match "\tT$" } | foreach { ($_ -split "`t")[0] } > ids_pos1.txt

# Position 2: 23349 (expect G)
seqkit subseq -r 23349:23349 $input | seqkit fx2tab | where { $_ -match "\tG$" } | foreach { ($_ -split "`t")[0] } > ids_pos2.txt

# Position 3: 28878 (expect T)
seqkit subseq -r 28878:28878 $input | seqkit fx2tab | where { $_ -match "\tT$" } | foreach { ($_ -split "`t")[0] } > ids_pos3.txt

# Position 4: 29591 (expect T)
seqkit subseq -r 29591:29591 $input | seqkit fx2tab | where { $_ -match "\tT$" } | foreach { ($_ -split "`t")[0] } > ids_pos4.txt

seqkit common ids_pos1.txt ids_pos2.txt ids_pos3.txt ids_pos4.txt > B.1.177_desired_ids.txt
seqkit grep -f B.1.177_desired_ids.txt $input > $output
