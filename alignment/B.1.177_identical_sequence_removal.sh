# Make output folder
New-Item -ItemType Directory E:\B117_dedup_by_sequence -Force | Out-Null

$in     = "E:\key_site_seqkit_filtered_B.1.177.dedup.fasta"
$outDir = "E:\B117_dedup_by_sequence"

# Output files
$outSeq = "$outDir\B.1.117_unique_bases_sequence.fasta"       # FASTA after removing identical sequences
$kept   = "$outDir\B.1.117.kept_unique_ids.txt"       # IDs that were kept
$removed= "$outDir\B.1.117.removed_identical_ids.txt"    # IDs that were removed

# 1) Remove identical sequences, write unique FASTA
seqkit rmdup -s -i -o $outSeq $in

# 2) Get kept IDs (headers from unique FASTA)
seqkit seq -n $outSeq > $kept

# 3) Get all IDs from input
seqkit seq -n $in > "$outDir\all_ids.txt"

# 4) Find removed IDs (all - kept)
Compare-Object (Get-Content "$outDir\all_ids.txt") (Get-Content $kept) -PassThru |
    Where-Object {$_ -ne ""} > $removed
