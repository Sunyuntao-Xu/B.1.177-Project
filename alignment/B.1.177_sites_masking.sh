# Begining of alignment start from 55's position of Wuhan-Hu-1

# T153G → Wuhan 153, align 99
seqkit subseq -r 99:99 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_T153G_check.tsv
Import-Csv E:/pos_T153G_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name

# G1149T → 1095
seqkit subseq -r 1095:1095 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_G1149T_check.tsv
Import-Csv E:/pos_G1149T_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name

# G2198A → 2144
seqkit subseq -r 2144:2144 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_G2198A_check.tsv
Import-Csv E:/pos_G2198A_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name

# G3145T → 3091
seqkit subseq -r 3091:3091 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_G3145T_check.tsv
Import-Csv E:/pos_G3145T_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name

# G3564T → 3510
seqkit subseq -r 3510:3510 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_G3564T_check.tsv
Import-Csv E:/pos_G3564T_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name

# A3778G → 3724
seqkit subseq -r 3724:3724 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_A3778G_check.tsv
Import-Csv E:/pos_A3778G_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name

# A4050C → 3996
seqkit subseq -r 3996:3996 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_A4050C_check.tsv
Import-Csv E:/pos_A4050C_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name

# C6255T → 6201
seqkit subseq -r 6201:6201 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_C6255T_check.tsv
Import-Csv E:/pos_C6255T_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name

# G5629 → 5575
seqkit subseq -r 5575:5575 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_5629G_check.tsv
Import-Csv E:/pos_5629G_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name

# A6851 → 6797
seqkit subseq -r 6797:6797 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_6851A_check.tsv
Import-Csv E:/pos_6851A_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name

# G7328 → 7274
seqkit subseq -r 7274:7274 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_7328G_check.tsv
Import-Csv E:/pos_7328G_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name

# T8022G → 7968
seqkit subseq -r 7968:7968 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_T8022G_check.tsv
Import-Csv E:/pos_T8022G_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name

# G8790T → 8736
seqkit subseq -r 8736:8736 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_G8790T_check.tsv
Import-Csv E:/pos_G8790T_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name

# T13402G → 13348
seqkit subseq -r 13348:13348 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_T13402G_check.tsv
Import-Csv E:/pos_T13402G_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name

# A13947T → 13893
seqkit subseq -r 13893:13893 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_A13947T_check.tsv
Import-Csv E:/pos_A13947T_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name

# C22802G → 22748
seqkit subseq -r 22748:22748 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_C22802G_check.tsv
Import-Csv E:/pos_C22802G_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name

# A24389C → 24335
seqkit subseq -r 24335:24335 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_A24389C_check.tsv
Import-Csv E:/pos_A24389C_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name

# G24390C → 24336
seqkit subseq -r 24336:24336 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_G24390C_check.tsv
Import-Csv E:/pos_G24390C_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name

# G24933T → 24879
seqkit subseq -r 24879:24879 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_G24933T_check.tsv
Import-Csv E:/pos_G24933T_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name

# A28095 → 28041
seqkit subseq -r 28041:28041 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_A28095_check.tsv
Import-Csv E:/pos_A28095_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name

# C29362 → 29308
seqkit subseq -r 29308:29308 E:/msaCodon_0201_B.1.177.final_trimmed.fasta | seqkit fx2tab > E:/pos_C29362_check.tsv
Import-Csv E:/pos_C29362_check.tsv -Delimiter "`t" -Header "ID", "Base" | Group-Object Base | Select-Object Count, Name
