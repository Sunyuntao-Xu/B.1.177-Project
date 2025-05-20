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
