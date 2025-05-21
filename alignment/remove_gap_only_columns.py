from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

# Step A: Define file paths
# Use raw strings to avoid issues with backslashes on Windows
input_file = r"E:\B.1.177_extracted.fasta"
output_file = r"E:\msaCodon_0201_B.1.177.cleaned.fasta"
print("Step A: File paths set.")

# Step B: Load the FASTA alignment
# Reads the multiple sequence alignment from the input file
alignment = AlignIO.read(input_file, "fasta")
print("Step B: Alignment loaded successfully.")

# Step C: Identify non-gap-only columns
# Iterates through each column and keeps only those that are not entirely gaps
filtered = []
for i in range(alignment.get_alignment_length()):
    column = alignment[:, i]
    if not all(c == '-' for c in column):
        filtered.append(i)
print(f"Step C: Found {len(filtered)} non-gap columns (of {alignment.get_alignment_length()} total).")

# Step D: Build new alignment without gap-only columns
# Constructs a new MultipleSeqAlignment object with only valid columns
new_alignment = MultipleSeqAlignment([
    record[:0] + ''.join(record.seq[i] for i in filtered)
    for record in alignment
])
print("Step D: New alignment created (gap-only columns removed).")

# Step E: Save cleaned alignment to output file
AlignIO.write(new_alignment, output_file, "fasta")
print(f"Step E: Cleaned alignment saved to: {output_file}")
