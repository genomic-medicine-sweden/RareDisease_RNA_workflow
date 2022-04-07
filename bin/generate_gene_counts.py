#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path
from collections import OrderedDict
import re


translator = {
    "unstranded": 1,
    "forward": 2,
    "reverse": 3,
}


def read_star_gene_cnts(sample: str, star: Path, strandedness: str) -> dict:
    """Read gene count file(s) from STAR output."""
    sample_ids = {}
    gene_ids = {}
    with open(star) as in_tab:
        for line in in_tab:
            if not line.startswith("N_"):
                gene_id = line.split()[0]
                strand = translator[strandedness]
                counts = line.split()[strand]
                gene_ids[gene_id] = int(counts)
    gene_ids = OrderedDict(sorted(gene_ids.items()))
    sample_ids[sample] = gene_ids
    return sample_ids


def transform_to_table(gene_ids_dict, outfile):
    """Transform in dictionary into tsv friendly."""
    rows = {}

    one_sample = next(iter(gene_ids_dict))
    gene_list = list(gene_ids_dict[one_sample].keys())

    for i, gene in enumerate(gene_list):
        sample_info = {}
        gene_id = {}
        for sample in gene_ids_dict:
            gene_id["GeneID"] = gene
            sample_info[sample] = gene_ids_dict[sample][gene]
            sample_info.update(gene_id)
        rows[i] = sample_info

    with open(outfile, "w", newline="") as tsvfile:
        sample_names = ["GeneID"] + [sample for sample in gene_ids_dict.keys()]
        writer = csv.DictWriter(tsvfile, fieldnames=sample_names, delimiter="\t")

        writer.writeheader()
        for keys in rows:
            writer.writerow(rows[keys])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.MetavarTypeHelpFormatter,
        description="""Generate collated gene counts from each STAR output.""",
    )
    parser.add_argument(
        "--star", type=str, nargs="+", help="*ReadsPerGene.out.tab from STAR"
    )
    parser.add_argument(
        "--sample", type=str, nargs="+", help="corresponding sample name"
    )
    parser.add_argument("--strandedness", type=str, help="strandedness of RNA")
    parser.add_argument("--output", type=str, help="output tsv file name")

    args = parser.parse_args()
    master_dict = {}
    for index, sample_id in enumerate(args.sample):
        sample_id = re.sub(r"[\[\],]", "", sample_id)
        master_dict.update(
            read_star_gene_cnts(
                sample=sample_id, star=args.star[index], strandedness=args.strandedness
            )
        )

    transform_to_table(master_dict, args.output)
