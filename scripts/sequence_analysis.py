
#!/usr/bin/env python3

"""
Project 1: Comparative mitochondrial COI gene analysis

Analyses real NCBI FASTA sequences:
- Sequence length
- GC content
- Nucleotide frequencies
- Sequence validation
- Pairwise sequence identity

Outputs:
- sequence_statistics.csv
- pairwise_identity.csv
- sequence_composition.png
- length_gc_content.png
"""

import csv
from pathlib import Path
from collections import Counter
from itertools import combinations

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from Bio import SeqIO
from Bio.Align import PairwiseAligner


# --------------------------------------------------
# 1. Project paths
# --------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent.parent

INPUT_FILE = BASE_DIR / "data/raw/all_species_COI.fasta"
RESULTS_DIR = BASE_DIR / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
PROCESSED_DIR = BASE_DIR / "data/processed"

for folder in [RESULTS_DIR, FIGURES_DIR, PROCESSED_DIR]:
    folder.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------
# 2. Load real DNA sequences
# --------------------------------------------------

if not INPUT_FILE.exists():
    raise FileNotFoundError(
        f"FASTA file not found: {INPUT_FILE}\n"
        "Run sequence_retrieval.py first."
    )

records = list(SeqIO.parse(INPUT_FILE, "fasta"))

if len(records) < 2:
    raise ValueError("At least two sequences are required.")

print(f"\nSequences loaded: {len(records)}")


# --------------------------------------------------
# 3. Sequence quality and composition analysis
# --------------------------------------------------

valid_bases = set("ACGTN")

statistics = []
sequences = {}

for record in records:

    sequence = str(record.seq).upper()

    if not sequence:
        print(f"WARNING: Empty sequence: {record.id}")
        continue

    invalid_bases = set(sequence) - valid_bases

    if invalid_bases:
        print(
            f"WARNING: {record.id} contains unsupported "
            f"characters: {invalid_bases}"
        )
        continue

    counts = Counter(sequence)

    length = len(sequence)
    gc_count = counts["G"] + counts["C"]


    valid_base_count = (
    counts["A"] + counts["T"] +
    counts["G"] + counts["C"]
    )
    gc_content = (
        gc_count / valid_base_count
    ) * 100 if valid_base_count else 0

    species = record.id.split("|")[0]

    sequences[species] = sequence

    statistics.append({
        "species": species,
        "accession": (
            record.id.split("|", 1)[1]
            if "|" in record.id else record.id
        ),
        "length_bp": length,
        "valid_base_count": valid_base_count,
        "A_count": counts["A"],
        "T_count": counts["T"],
        "G_count": counts["G"],
        "C_count": counts["C"],
        "N_count": counts["N"],
        "GC_content_percent": round(gc_content, 3),
        "ambiguous_base_percent": round(
            (counts["N"] / length) * 100, 3
        )
    })

if len(sequences) < 2:
    raise ValueError("Not enough valid sequences for analysis.")

df = pd.DataFrame(statistics)

output_stats = PROCESSED_DIR / "sequence_statistics.csv"
df.to_csv(output_stats, index=False)

print("\nSEQUENCE STATISTICS")
print(df.to_string(index=False))


# --------------------------------------------------
# 4. Pairwise sequence identity
# --------------------------------------------------

aligner = PairwiseAligner()

aligner.mode = "global"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -5
aligner.extend_gap_score = -0.5

pairwise_results = []

for (species_a, seq_a), (species_b, seq_b) in combinations(
    sequences.items(), 2
):

    alignment = aligner.align(seq_a, seq_b)[0]
    indices = alignment.indices

    comparable_sites = 0
    identical_sites = 0
    gap_sites = 0

    for index_a, index_b in zip(indices[0], indices[1]):

        if index_a == -1 or index_b == -1:
            gap_sites += 1
            continue

        base_a = seq_a[index_a]
        base_b = seq_b[index_b]

        if base_a in "ACGT" and base_b in "ACGT":
            comparable_sites += 1

            if base_a == base_b:
                identical_sites += 1

    identity = (
        identical_sites / comparable_sites * 100
        if comparable_sites else 0
    )

    pairwise_results.append({
        "species_1": species_a,
        "species_2": species_b,
        "identical_sites": identical_sites,
        "comparable_sites": comparable_sites,
        "gap_sites": gap_sites,
        "identity_percent": round(identity, 3)
    })

    print(
        f"{species_a} vs {species_b}: "
        f"{identity:.2f}% identity"
    )

pairwise_df = pd.DataFrame(pairwise_results)

pairwise_file = PROCESSED_DIR / "pairwise_identity.csv"
pairwise_df.to_csv(pairwise_file, index=False)

# --------------------------------------------------
# 5. Nucleotide composition plot
# --------------------------------------------------

base_columns = ["A_count", "T_count", "G_count", "C_count"]

plot_df = df.set_index("species")[base_columns]

ax = plot_df.plot(
    kind="bar",
    figsize=(10, 6)
)

ax.set_title("Nucleotide Composition of COI Sequences")
ax.set_xlabel("Species")
ax.set_ylabel("Nucleotide Count")
ax.legend(title="Nucleotide")
plt.xticks(rotation=30, ha="right")
plt.tight_layout()

plt.savefig(
    FIGURES_DIR / "sequence_composition.png",
    dpi=300
)

plt.close()


# --------------------------------------------------
# 6. Sequence length and GC content plot
# --------------------------------------------------

fig, ax1 = plt.subplots(figsize=(10, 6))

x = np.arange(len(df))

ax1.bar(
    x,
    df["length_bp"],
    alpha=0.7,
    label="Sequence length (bp)"
)

ax1.set_ylabel("Sequence length (bp)")
ax1.set_xticks(x)
ax1.set_xticklabels(df["species"], rotation=30, ha="right")

ax2 = ax1.twinx()

ax2.plot(
    x,
    df["GC_content_percent"],
    marker="o",
    linestyle="-",
    label="GC content (%)"
)

ax2.set_ylabel("GC content (%)")

ax1.set_title("COI Sequence Length and GC Content")

fig.tight_layout()

plt.savefig(
    FIGURES_DIR / "length_gc_content.png",
    dpi=300
)

plt.close()


# --------------------------------------------------
# 7. Final report
# --------------------------------------------------

print("\n" + "=" * 55)
print("COMPARATIVE COI ANALYSIS COMPLETED")
print("=" * 55)

print(f"Valid sequences analysed: {len(sequences)}")
print(f"Statistics: {output_stats}")
print(f"Pairwise identity: {pairwise_file}")

print(f"Figures: {FIGURES_DIR}")

print("\nAnalysis outputs:")
print("- sequence_statistics.csv")
print("- pairwise_identity.csv")
print("- sequence_composition.png")
print("- length_gc_content.png")
