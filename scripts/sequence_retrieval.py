#!/usr/bin/env python3

"""
Project 1: Comparative COI Gene Sequence Analysis

Retrieve mitochondrial COI coding sequences from NCBI
for five vertebrate species.

Data source: NCBI Nucleotide / GenBank
"""

import os
import csv
import time
from pathlib import Path

from Bio import Entrez, SeqIO


# --------------------------------------------------
# 1. Configuration
# --------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent.parent

RAW_DIR = BASE_DIR / "data" / "raw"
PROCESSED_DIR = BASE_DIR / "data" / "processed"

RAW_DIR.mkdir(parents=True, exist_ok=True)
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

SPECIES = {
    "human": "Homo sapiens",
    "mouse": "Mus musculus",
    "dog": "Canis lupus familiaris",
    "cattle": "Bos taurus",
    "chicken": "Gallus gallus",
}

email = os.environ.get("NCBI_EMAIL")

if not email:
    raise ValueError(
        "NCBI_EMAIL is missing. Set your email before running this script."
    )

Entrez.email = email
Entrez.tool = "Comparative_COI_Analysis"


# --------------------------------------------------
# 2. Find annotated COI coding sequences
# --------------------------------------------------

def find_coi_sequence(record):

    candidates = []

    for feature in record.features:

        if feature.type != "CDS":
            continue

        qualifiers = feature.qualifiers

        gene_names = [
            value.lower()
            for value in qualifiers.get("gene", [])
        ]

        products = [
            value.lower()
            for value in qualifiers.get("product", [])
        ]

        is_coi = (
            any(name in ["coi", "cox1", "cox-1"] for name in gene_names)
            or any(
                "cytochrome c oxidase subunit i" in product
                for product in products
            )
        )

        if not is_coi:
            continue

        if "pseudo" in qualifiers or "pseudogene" in qualifiers:
            continue

        sequence = feature.extract(record.seq)

        if 1000 <= len(sequence) <= 1800:
            candidates.append((len(sequence), sequence))

    if not candidates:
        return None

    candidates.sort(key=lambda item: item[0], reverse=True)

    return candidates[0][1]


# --------------------------------------------------
# 3. Retrieve sequences from NCBI
# --------------------------------------------------

def retrieve_species_sequence(label, scientific_name):

    print(f"\nSearching NCBI for {scientific_name}...")

    query = (
        f'"{scientific_name}"[Organism] AND '
        '(COI[Gene] OR COX1[Gene]) AND '
        'mitochondrion[filter]'
    )

    try:
        search_handle = Entrez.esearch(
            db="nuccore",
            term=query,
            retmax=20,
            sort="relevance"
        )

        search_results = Entrez.read(search_handle)
        search_handle.close()

        record_ids = search_results["IdList"]

        if not record_ids:
            print(f"No records found for {scientific_name}")
            return None

        fetch_handle = Entrez.efetch(
            db="nuccore",
            id=",".join(record_ids),
            rettype="gb",
            retmode="text"
        )

        records = list(SeqIO.parse(fetch_handle, "genbank"))
        fetch_handle.close()

        for record in records:

            sequence = find_coi_sequence(record)

            if sequence is None:
                continue

            record_id = record.id

            output_record = record[:]
            output_record.seq = sequence
            output_record.id = f"{label}|{record_id}"
            output_record.name = label
            output_record.description = (
                f"{scientific_name} mitochondrial COI"
            )

            output_file = RAW_DIR / f"{label}_COI.fasta"

            SeqIO.write(output_record, output_file, "fasta")

            print(f"Selected accession: {record_id}")
            print(f"Sequence length: {len(sequence)} bp")
            print(f"Saved: {output_file}")

            return {
                "species": scientific_name,
                "label": label,
                "accession": record_id,
                "length_bp": len(sequence),
                "sequence_file": str(output_file.relative_to(BASE_DIR))
            }

        print(f"No suitable annotated COI CDS found for {scientific_name}")
        return None

    except Exception as error:
        print(f"Error retrieving {scientific_name}: {error}")
        return None


# --------------------------------------------------
# 4. Run retrieval for all species
# --------------------------------------------------

def main():

    metadata = []

    for label, scientific_name in SPECIES.items():

        result = retrieve_species_sequence(label, scientific_name)

        if result:
            metadata.append(result)

        # Keep requests within NCBI's usage guidance.
        time.sleep(0.5)

    if not metadata:
        print("\nNo sequences were retrieved.")
        return

    # Save metadata
    metadata_file = PROCESSED_DIR / "sequence_metadata.csv"

    with open(metadata_file, "w", newline="") as csvfile:

        writer = csv.DictWriter(
            csvfile,
            fieldnames=[
                "species",
                "label",
                "accession",
                "length_bp",
                "sequence_file"
            ]
        )

        writer.writeheader()
        writer.writerows(metadata)

    # Create combined FASTA
    fasta_files = [
        RAW_DIR / f"{label}_COI.fasta"
        for label in SPECIES
        if (RAW_DIR / f"{label}_COI.fasta").exists()
    ]

    combined_file = RAW_DIR / "all_species_COI.fasta"

    with open(combined_file, "w") as output:

        for fasta_file in fasta_files:
            with open(fasta_file) as input_file:
                output.write(input_file.read())

    print("\n----------------------------------")
    print("SEQUENCE RETRIEVAL COMPLETE")
    print("----------------------------------")
    print(f"Sequences retrieved: {len(metadata)}/{len(SPECIES)}")
    print(f"Metadata: {metadata_file}")
    print(f"Combined FASTA: {combined_file}")


if __name__ == "__main__":
    main()
