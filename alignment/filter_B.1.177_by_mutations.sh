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

$allBad = $bad1 + $bad2 + $bad3 + $bad4 | Sort-Object -Unique
$allBad > "E:/B.1.177_exclude_ids.txt"
