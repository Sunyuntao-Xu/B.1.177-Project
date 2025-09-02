# ================================
# Deduplicate B.1.177 by sequence
# ================================

# --- Setup ---
New-Item -ItemType Directory E:\B117_dedup_by_sequence -Force | Out-Null

$in     = "E:\key_site_seqkit_filtered_B.1.177.dedup.fasta"
$outDir = "E:\B117_dedup_by_sequence"

# Output files
$outSeq = "$outDir\B117.unique.fasta"                # final FASTA (unique by bases)
$dupFa  = "$outDir\B117.duplicates.fa.gz"            # removed duplicate records
$dupMap = "$outDir\B117.identical.detail.txt"        # duplicate_id<TAB>kept_id
$log    = "$outDir\B117.sequence_dedup.log.txt"

# --- 1) Input stats ---
seqkit stats $in | Tee-Object -FilePath $log
"Seq count (before): " + ((Get-Content $in | Select-String '^>').Count) |
    Tee-Object -FilePath $log -Append

# --- 2) Remove identical sequences (by bases) ---
# -s : dedupe by sequence content
# -j : threads
# -D : mapping duplicate_id -> kept_id
# -d : save duplicates
# seqkit seq -w 0 : disable line wrapping
seqkit rmdup -s -j 8 -D $dupMap -d $dupFa $in |
    seqkit seq -w 0 |
    Out-File -Encoding ascii $outSeq

# --- 3) Output stats ---
"Seq count (after):  "  + ((Get-Content $outSeq | Select-String '^>').Count) |
    Tee-Object -FilePath $log -Append
"Duplicates mapped:  "  + ((Get-Content $dupMap).Count) |
    Tee-Object -FilePath $log -Append
"`nFirst few mappings (duplicate_id -> kept_id):" |
    Tee-Object -FilePath $log -Append
Get-Content $dupMap -TotalCount 10 | Tee-Object -FilePath $log -Append

# --- 4) Sanity check (print first header & length) ---
"First record in unique FASTA:" | Tee-Object -FilePath $log -Append
seqkit head -n 1 $outSeq | Tee-Object -FilePath $log -Append
