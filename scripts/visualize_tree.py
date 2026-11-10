from Bio import Phylo
import matplotlib.pyplot as plt
from pathlib import Path

# Input and output paths
tree_file = Path("results/phylogenetic_trees/COI_tree.nwk")
output_file = Path("results/figures/COI_phylogenetic_tree.png")

# Load tree
tree = Phylo.read(tree_file, "newick")
tree.rooted = False

# Draw tree
fig = plt.figure(figsize=(12, 7))

axes = fig.add_subplot(1, 1, 1)

Phylo.draw(
    tree,
    axes=axes,
    do_show=False,
    show_confidence=True
)

axes.set_title(
    "Phylogenetic Analysis of COI Sequences",
    fontsize=14
)

axes.set_xlabel("Evolutionary distance (substitutions per site)")

plt.tight_layout()
plt.savefig(output_file, dpi=300, bbox_inches="tight")
plt.close()

print(f"Tree visualization saved: {output_file}")
