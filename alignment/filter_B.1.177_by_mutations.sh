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


# Fix the GISAID ID list by replacing spaces with underscores
Get-Content "E:\GISAID_selected_sequences.txt" |
    ForEach-Object { $_.Trim() -replace " ", "_" } |
    Where-Object { $_ -ne "" } |
    Set-Content "E:\GISAID_selected_sequences_nospace.txt"


Get-Content "E:\B.1.177_exclude_ids.txt" |
    ForEach-Object { $_.Trim() -replace " ", "_" } |
    Where-Object { $_ -ne "" } |
    Set-Content "E:\B.1.177_exclude_ids_nospace.txt"

$all = Get-Content "E:\GISAID_selected_sequences_nospace.txt"
$bad = Get-Content "E:\B.1.177_exclude_ids_nospace.txt"

# Continue as before
$badSet = [System.Collections.Generic.HashSet[string]]::new()
$bad | ForEach-Object { $null = $badSet.Add($_) }

$wanted = $all | Where-Object { -not $badSet.Contains($_) }
$wanted | Sort-Object -Unique > "E:\B.1.177_desired_ids.txt"


# Set file paths
$fastaFile = "E:\masked_output.fasta"                    # ✅ Use the masked FASTA
$nameListFile = "E:\B.1.177_desired_ids.txt"             # ✅ List of wanted sequence IDs (with _ not spaces)
$outputFile = "E:\key_site_filtered_B.1.177.fasta"       # ✅ Output file

# Load target headers into hashset
$targetNames = Get-Content $nameListFile | ForEach-Object { $_.Trim() } | Where-Object { $_ -ne "" }
$targetSet = @{}
foreach ($name in $targetNames) {
    $targetSet[$name] = $true
}

# Read masked FASTA and extract matching sequences
$writeBlock = $false
Get-Content $fastaFile | ForEach-Object {
    if ($_ -like ">*") {
        $header = $_.Substring(1).Trim()
        if ($targetSet.ContainsKey($header)) {
            $writeBlock = $true
            Add-Content -Path $outputFile -Value $_
        } else {
            $writeBlock = $false
        }
    } elseif ($writeBlock) {
        Add-Content -Path $outputFile -Value $_
    }
}
