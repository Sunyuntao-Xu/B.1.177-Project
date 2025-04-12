# Split the cleaned B.1.177 alignment into 20 chunks to reduce memory load for MAFFT
# Each chunk will be aligned independently to the reference to avoid overloading RAM
seqkit split -p 20 E:\msaCodon_0201_B.1.177.final_trimmed.fasta -O E:\split_input\

# Open WSL to run MAFFT in a Linux-like environment
wsl

# Navigate to the directory where split FASTA files are located (in WSL path)
cd /mnt/e/split_input

# For each split FASTA chunk:
# - Use MAFFT with --add and --keeplength to add it to the reference alignment
# - --add ensures new sequences align to existing reference
# - --keeplength preserves the length and structure of the reference alignment
# - After each addition, overwrite the reference with the updated alignment
# This incremental approach allows large-scale alignments to complete with limited memory
for file in *.fasta; do
    echo "Adding $file to reference..."
    mafft --add "$file" \
          --keeplength "/mnt/e/Reference Alignments/merged_reference_sequences_aligned_masked.fasta" \
          > "/mnt/e/temp_output.fasta"

    # Update the reference alignment with newly added sequences
    cp /mnt/e/temp_output.fasta "/mnt/e/Reference Alignments/merged_reference_sequences_aligned_masked.fasta"
done
