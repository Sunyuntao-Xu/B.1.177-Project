New-Item -ItemType Directory E:\B117_dedup_by_sequence -Force | Out-Null

$in     = "E:\key_site_seqkit_filtered_B.1.177.dedup.fasta"
$outDir = "E:\B117_dedup_by_sequence"

# Core outputs
$outSeq = "$outDir\B117.unique.fasta"                # unique by sequence (exact bases)
$dupFa  = "$outDir\B117.duplicates.fa.gz"            # removed duplicate records
$dupMap = "$outDir\B117.identical.detail.txt"        # duplicate_id<TAB>kept_id mapping
$log    = "$outDir\B117.sequence_dedup.log.txt"


seqkit stats $in | Tee-Object -FilePath $log
"Seq count (before): " + ((Get-Content $in | Select-String '^>').Count) | Tee-Object -FilePath $log -Append


# -s : dedupe by sequence (not name)
# -D : write mapping duplicate_id<TAB>kept_id
# -d : save removed duplicates as FASTA (gz)
# -j : threads
seqkit rmdup -s -j 8 -D $dupMap -d $dupFa $in > $outSeq

"Seq count (after):  "  + ((Get-Content $outSeq | Select-String '^>').Count) | Tee-Object -FilePath $log -Append
"Duplicates mapped:  "  + ((Get-Content $dupMap).Count)                    | Tee-Object -FilePath $log -Append
"`nFirst few mappings (duplicate_id -> kept_id):"                           | Tee-Object -FilePath $log -Append
Get-Content $dupMap -TotalCount 10 | Tee-Object -FilePath $log -Append
