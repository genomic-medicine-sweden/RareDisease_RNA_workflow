#!/usr/bin/env python3

import argparse
import os
import csv

d = {
    "RNA_ID": "",
    "RNA_BAM_FILE": "",
    "DNA_VCF_FILE": "",
    "DNA_ID": "",
    "DROP_GROUP": "blood",
    "PAIRED_END": "",
    "COUNT_MODE": "",
    "COUNT_OVERLAPS": "",
    "STRAND": "",
    "HPO_TERMS": "",
    "GENE_COUNTS_FILE": "",
    "GENE_ANNOTATION": "",
    "GENOME": "",
}


def populate_dictionary(rna_id, cnts):

    master_dict = {}
    for i, rna in enumerate(rna_id):
        temp_d = d.copy()
        temp_d["RNA_ID"] = rna

        if cnts:
            temp_d["GENE_COUNTS_FILE"] = cnts
            temp_d["GENE_ANNOTATION"] = "hg38"

        master_dict[i] = temp_d

    return master_dict


def write_table(dict, outfile):
    with open(outfile, "w", newline="") as tsvfile:
        header = dict[0].keys()
        writer = csv.DictWriter(tsvfile, fieldnames=header, delimiter="\t")

        writer.writeheader()
        for index in dict:
            writer.writerow(dict[index])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.MetavarTypeHelpFormatter,
        description="""Generate sample annotation file for DROP.""",
    )

    parser.add_argument(
        "--samples", type=str, help="A list of RNA sample names delimited by space."
    )
    parser.add_argument(
        "--counts",
        type=str,
        help="A tsv file of gene counts for all processed samples.",
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Path to save to",
    )

    args = parser.parse_args()
    sample_list = [str(item) for item in args.samples.split()]

    write_table(populate_dictionary(sample_list, args.counts), args.output)
