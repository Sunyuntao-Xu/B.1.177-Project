# Position 1: 22227 (C22227T → expect T)
seqkit subseq -r 22173:22173 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_22227_check.tsv
Import-Csv E:/pos_22227_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name

# Position 2: 23403 (A23403G → expect G)
seqkit subseq -r 23349:23349 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_23403_check.tsv
Import-Csv E:/pos_23403_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name

# Position 3: 28932 (C28932T → expect T)
seqkit subseq -r 28878:28878 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_28932_check.tsv
Import-Csv E:/pos_28932_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name

# Position 4: 29654 (G29654T → expect T or G depending on variant)
seqkit subseq -r 29600:29600 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_29654_check.tsv
Import-Csv E:/pos_29654_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name
