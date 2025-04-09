from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

print("a")

# Use raw string to handle backslashes in Windows paths
input_file = r"E:\msaCodon_0201_B.1.177.fasta"
output_file = r"E:\msaCodon_0201_B.1.177.cleaned.fasta"

print("b")

# Load the alignment
alignment = AlignIO.read(input_file, "fasta")

print("c")

# Identify columns that are not gap-only
filtered = []
for i in range(alignment.get_alignment_length()):
    column = alignment[:, i]
    if not all(c == '-' for c in column):
        filtered.append(i)

# Create a new alignment without gap-only columns
new_alignment = MultipleSeqAlignment([
    record[:0] + ''.join(record.seq[i] for i in filtered)
    for record in alignment
])

print("d")

# Write the cleaned alignment to a new file
AlignIO.write(new_alignment, output_file, "fasta")

print("e")

print("Gap-only columns removed. Cleaned file saved to:")
print(output_file)

print("f")
