#!/usr/bin/env python
# Python 3.7 
import argparse
import pandas as pd

# Global variables
fasta_lines = []
gtf_lines = []

# Taken from crop-seq
fasta_header_template = ">{chrom} dna:chromosome chromosome:GRCh38:{chrom}:1:{length}:1 REF"

gtf_template = """{chrom}\thavana\tgene\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; gene_name "{id}_gene"; gene_source "ensembl_havana"; gene_biotype "lincRNA";
{chrom}\thavana\ttranscript\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; transcript_id "{id}_transcript"; gene_name "{id}_gene"; gene_source "ensembl_havana"; gene_biotype "lincRNA"; transcript_name "{id}_transcript"; transcript_source "havana";
{chrom}\thavana\texon\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; transcript_id "{id}_transcript"; exon_number "1"; gene_name "{id}_gene"; gene_source "ensembl_havana"; gene_biotype "lincRNA"; transcript_name "{id}_transcript"; transcript_source "havana"; exon_id "{id}_exon";
"""

# Build sequence


def buildRef(row):
    ID = row['ID']
    ref_seq = row["Sequence"]
    header = fasta_header_template.format(chrom=ID, length=len(ref_seq))
    fasta_lines.append(header)
    fasta_lines.append(ref_seq)
    gtf_lines.append(gtf_template.format(chrom=ID, id=ID, length=len(ref_seq)))


# Parsers
def loadFiles(args):
    cropseq_df = pd.read_csv(args.cropseq)
    return(cropseq_df)

# Argument parser
def parseArgs():
    parser = argparse.ArgumentParser(prog = "PrepRef")
    parser.add_argument("-c", "--cropseq", type = str, help = "CROP-seq file")
    parser.add_argument("-o", "--output", type = str, help = "Basename of output files")
    args = parser.parse_args()
    return(args)

if __name__ ==  "__main__":
    # Parse arguments
    args = parseArgs()

    # Parse spreadsheets
    cropseq_df = loadFiles(args)

    cropseq_df.apply(buildRef, axis = 1)

    # write to file
    output_fasta = args.output + ".fa"
    output_gtf = args.output + ".gtf"
    
    with open(output_fasta, "w") as fasta_handle:
        fasta_handle.writelines("\n".join(fasta_lines))
    with open(output_gtf, "w") as gtf_handle:
        gtf_handle.writelines(gtf_lines)




