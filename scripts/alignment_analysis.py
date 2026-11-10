from Bio import AlignIO
import pandas as pd
from pathlib import Path

# Project paths
ALIGNMENT_FILE = Path("results/alignments/COI_aligned.fasta")
OUTPUT_DIR = Path("data/processed")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Load alignment
alignment = AlignIO.read(ALIGNMENT_FILE, "fasta")

total_sites = alignment.get_alignment_length()
total_sequences = len(alignment)

conserved_sites = 0
variable_sites = 0
parsimony_informative_sites = 0
ambiguous_sites = 0

# Analyze each alignment column
for column in range(total_sites):

    bases = [
        record.seq[column].upper()
        for record in alignment
    ]

    # Exclude gaps and ambiguous nucleotides
    valid_bases = [
        base for base in bases
        if base in "ACGT"
    ]

    # Count columns containing gaps or ambiguous bases
    if len(valid_bases) != total_sequences:
        ambiguous_sites += 1
        continue

    counts = pd.Series(valid_bases).value_counts()

    if len(counts) == 1:
        conserved_sites += 1

    else:
        variable_sites += 1

        if (counts >= 2).sum() >= 2:
            parsimony_informative_sites += 1

# Results
results = {
    "total_alignment_sites": total_sites,
    "conserved_sites": conserved_sites,
    "variable_sites": variable_sites,
    "parsimony_informative_sites": parsimony_informative_sites,
    "ambiguous_or_gap_sites": ambiguous_sites
}

# Save results
results_df = pd.DataFrame(
    list(results.items()),
    columns=["Parameter", "Value"]
)

results_df.to_csv(
    OUTPUT_DIR / "alignment_statistics.csv",
    index=False
)

print("\nMULTIPLE SEQUENCE ALIGNMENT ANALYSIS")
print("-----------------------------------")

for parameter, value in results.items():
    print(f"{parameter}: {value}")

print("\nResults saved to data/processed/alignment_statistics.csv")
