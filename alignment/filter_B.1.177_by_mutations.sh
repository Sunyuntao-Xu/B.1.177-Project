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
$fastaFile = "E:\masked_output.fasta"
$nameListFile = "E:\B.1.177_desired_ids.txt"
$outputFile = "E:\key_site_filtered_B.1.177.fasta"

# Load desired sequence names into a HashSet
$targetSet = [System.Collections.Generic.HashSet[string]]::new()
Get-Content $nameListFile | ForEach-Object {
    $name = $_.Trim()
    if ($name -ne "") { $null = $targetSet.Add($name) }
}

# Open FASTA reader and writer
$reader = [System.IO.StreamReader]::new($fastaFile)
$writer = [System.IO.StreamWriter]::new($outputFile, $false, [System.Text.Encoding]::UTF8)
$writer.NewLine = "`n"  # Use Unix line endings for compact output

# Extract matching sequences
$writeBlock = $false

while (-not $reader.EndOfStream) {
    $line = $reader.ReadLine()

    if ($line.StartsWith(">")) {
        $header = $line.Substring(1).Trim()
        $writeBlock = $targetSet.Contains($header)
        if ($writeBlock) { $writer.WriteLine($line) }
    } elseif ($writeBlock) {
        $writer.WriteLine($line)
    }
}

# Close the streams
$reader.Close()
$writer.Close()

Write-Output "✅ FASTA extraction complete: $outputFile"

