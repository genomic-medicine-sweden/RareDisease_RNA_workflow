#!/usr/bin/env python3

import argparse
import csv
import os
from collections import OrderedDict


translator = {
    "unstranded": 1,
    "forward": 2,
    "reverse": 3,
}


def read_star_gene_cnts(sample, star, strandedness):
    """Read gene counts file from STAR output."""
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


def read_existing_gene_cnts(path_tsv, strandedness):
    """Read existing gene counts file."""
    sample_ids = {}
    gene_ids = {}
    with open(path_tsv) as in_tsv:
        for line in in_tsv:
            if line.startswith("GeneID"):
                sample = line.split()[-1]
            else:
                gene_id = line.split()[0]
                strand = translator[strandedness]
                counts = line.split()[strand]
                gene_ids[gene_id] = int(counts)
    gene_ids = OrderedDict(sorted(gene_ids.items()))
    sample_ids[sample] = gene_ids
    return sample_ids


def combine_dict(dict_1, dict_2):
    """Combine the existing file->dict to the incoming one."""
    master = {}
    master.update(dict_1)
    master.update(dict_2)
    return master


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


def pad_table(d1, d2):
    """Pad two different sized dictionaries."""
    # take one sample from each dictionary
    d1_sample = next(iter(d1))
    d2_sample = next(iter(d2))

    # get the genes to compare
    d1_list = set(d1[d1_sample].keys())
    d2_list = set(d2[d2_sample].keys())

    # aggregate genes from both lists and take unique
    master_list = set(list(d1_list) + list(d2_list))

    for gene in master_list - d2_list:
        d2[d2_sample][gene] = 0
        d2_list.add(gene)

    master_dict = combine_dict(d1, d2)
    return master_dict


def file_exists(path_tsv):
    """Check if file exists."""
    if os.path.exists(path_tsv):
        return path_tsv
    else:
        return False


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.MetavarTypeHelpFormatter,
        description="""Generate collated gene counts from each STAR output.""",
    )
    parser.add_argument("--star", type=str, help="*ReadsPerGene.out.tab from STAR")
    parser.add_argument("--sample", type=str, help="corresponding sample name")
    parser.add_argument("--strandedness", type=str, help="strandedness of RNA")
    parser.add_argument("--output", type=str, help="output tsv file name")

    args = parser.parse_args()
    file_exist = file_exists("external_geneCounts.tsv")
    in_dict = read_star_gene_cnts(args.sample, args.star, args.strandedness)
    if file_exist:
        transform_to_table(
            pad_table(in_dict, read_existing_gene_cnts(file_exist), args.strandedness),
            args.output,
        )
    else:
        transform_to_table(in_dict, args.output)
