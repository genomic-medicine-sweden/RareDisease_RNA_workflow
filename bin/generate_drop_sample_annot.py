#!/usr/bin/env python3

import argparse
import os
import csv


class SampleAnnotation:
    """SampleAnnotation class"""

    SAMPLE_ANNOTATION_COLUMNS = [
        "RNA_ID",
        "RNA_BAM_FILE",
        "DNA_VCF_FILE",
        "DNA_ID",
        "DROP_GROUP",
        "PAIRED_END",
        "COUNT_MODE",
        "COUNT_OVERLAPS",
        "STRAND",
        "HPO_TERMS",
        "GENE_COUNTS_FILE",
        "GENE_ANNOTATION",
        "GENOME",
    ]

    def __init__(self, cnts_file, out_file):
        """Create SampleAnnotation given the parameters"""
        self.cnts_file = cnts_file
        self.out_file = out_file

    def parse_header(self) -> list[str]:
        """Parse the first line of gene counts file"""
        header = []
        with open(self.cnts_file) as file_object:
            header = file_object.readline().split()

        del header[0]  # remove GeneID field
        return header

    def write_table(self):
        """Write the Sample Annotation tsv file"""
        with open(self.out_file, "w") as tsv_file:
            fieldnames = self.SAMPLE_ANNOTATION_COLUMNS
            sample_ids = self.parse_header()
            writer = csv.DictWriter(tsv_file, fieldnames=fieldnames, delimiter="\t")

            writer.writeheader()
            for id in sample_ids:
                sa_dict = {}.fromkeys(fieldnames, "")
                sa_dict["RNA_ID"] = id
                sa_dict["DROP_GROUP"] = "blood"
                sa_dict["GENE_COUNTS_FILE"] = os.path.basename(self.cnts_file)
                writer.writerow(sa_dict)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.MetavarTypeHelpFormatter,
        description="""Generate sample annotation file for DROP.""",
    )

    parser.add_argument(
        "--count_file",
        type=str,
        help="A tsv file of gene counts for all processed samples.",
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Path to save to",
    )

    args = parser.parse_args()

    SampleAnnotation(args.counts, args.output).write_table()
