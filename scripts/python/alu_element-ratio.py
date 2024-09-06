#!/usr/bin/env python3
#
# Script Name: Alu_element-ratio.py
#
# Author: Florian Janke
# Last updated: 2024/09/06
#
# Description
#   (1) Count the number of reads spanning (= non-aberrant) or overlapping recurrently protected regions (= aberrant)
#   (2) Returns results per GC content and fragment length
#
#---------------------------------------------------------------------------------------------------

#-----------
# P A C K A G E S
#-----------
import numpy as np
import os
import argparse
from multiprocessing import Pool
from functools import partial
import pandas as pd
import pysam
from Bio.SeqUtils import GC
from FrEIA_tools import CastDataTypes


#-----------
# P A R A M E T E R S
#-----------
def ParsingArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile",
                        dest="infile",
                        type=str,
                        required=True,
                        help="Path to the input bam file.")
    parser.add_argument("-o", "--outfile",
                        dest="outfile",
                        type=str,
                        required=True,
                        help="Path to output .tsv file.")
    parser.add_argument("-x", "--regions",
                        dest="regions",
                        type=str,
                        required=True,
                        help="Path to regions to be processed .bed file.")
    parser.add_argument("-y", "--type",
                        dest="type",
                        type=str,
                        required=True,
                        help="Type of regions to bed processed. Currently supported: Alu, CpG_islands, Whole_genome")
    parser.add_argument("-c", "--contigs",
                        dest="contigs",
                        type=str,
                        required=True,
                        help="List of contig names separated by a coma.")
    parser.add_argument("-r", "--reference",
                        dest="ref",
                        type=str,
                        required=True,
                        help="Path to genome reference .fasta file.")
    parser.add_argument("-t", "--threads",
                        dest="threads",
                        type=int,
                        required=True,
                        help="Number of threads.")
    return parser.parse_args()


#-----------
# F U N C T I O N S
#-----------
def threading_parts(args, chunk):
    # Load chromosomes to partition
    chr = args.contigs.split(",")
    # Load bam file
    bamfile = pysam.AlignmentFile(args.infile, "rb")
    for i in chr:
        # Extract chromosome length
        chrom_lengths = dict(zip(bamfile.references, bamfile.lengths))[i]
        # Determine number of parts to get the required <chunk> size
        part = int(round(chrom_lengths /chunk + 0.5, 0))
        # Define chromosome chunk
        chr_chunk = int(round(chrom_lengths /part, 0))
        # Get chromosomal coordinates of chunks
        for n in range(1, part+1):
            if n == 1:
                start = 1
            else:
                start = end +1
            end = start +chr_chunk -1
            if n == part:
                end = chrom_lengths
            if n == 1:
                coord = [i + '_' + str(start) + '_' + str(end)]
            else:
                coord = coord + [i + '_' + str(start) + '_' + str(end)]
        # Append chunks of chromosomes
        if i == chr[0]:
            chromosome_parts = coord
        else:
            chromosome_parts = chromosome_parts + coord
    bamfile.close()
    # Return <chromosome_parts>
    return chromosome_parts

def subset_rpr(args, chromosomes):
    # Create output folder
    chunk_dir = os.path.join(os.path.dirname(args.regions), args.type)
    os.makedirs(chunk_dir, exist_ok=True)  
    # Check whether all files were generated already
    for file in chromosomes:
        file_path = os.path.join(chunk_dir, f"{file}.npz")
        tmp = [os.path.exists(file_path)]
        if file == chromosomes[0]:
            statement = tmp
        else:
            statement = statement + tmp
    if all(statement):
        return print('.npz were already generated.')
    # Load regions
    regions = pd.read_csv(args.regions, sep="\t")
    for chunk in chromosomes:
        chunk = chunk.split('_')
        # Skip if file already existis
        if os.path.exists(os.path.join(chunk_dir, f"{chunk[0]}_{chunk[1]}_{chunk[2]}.npz")):
            continue
        # Subset
        tmp = regions[regions['chr'] == chunk[0]]
        if int(chunk[1]) -1000 < 0:
            start = int(chunk[1])
        else:
            start = int(chunk[1]) -1000
        end = int(chunk[2]) +1000
        tmp = tmp[(tmp['start'] >= start) & (tmp['end'] < end)]
        tmp = tmp[tmp.columns.drop('chr')]
        tmp = tmp.to_numpy()
        # Save
        np.savez(os.path.join(chunk_dir, f"{chunk[0]}_{chunk[1]}_{chunk[2]}.npz"), my_array=tmp)

def fetchData_tumorSpecific(args, chromosome):
    # Define chromosome, start and end position
    chrom = chromosome
    # Initiate output dictionary
    outDict = {"ID": [],
               "CGN": [],
               "NCG": []}
    # Load .sam file
    samfile = pysam.FastaFile(args.ref)
    # Conntect .bam file to pysam library
    bamfile = pysam.AlignmentFile(args.infile, "rb")
    regions = pd.read_csv(args.regions, sep="\t", header=None)
    regions = regions[regions.iloc[:, 0] == chrom]
    regions = regions.iloc[:, 1].values
    regions = regions -1
    for read in bamfile.fetch(chrom):
        # Skip read that are not paired, properly paired, have low MAPQ value and are hard-clipped
        if not read.is_paired or not read.is_proper_pair or read.mapping_quality < 10 or 'H' in read.cigarstring or read.is_reverse:
            continue
        # Extract 5' end motif
        motif = read.query_alignment_sequence[:3]
        # Check if the end motif equals CGN or NCG
        if read.pos in regions:
            if motif[:2].upper() == "CG":
                type = "CGN"
            else:
                continue
        elif read.pos in regions-1:
            if motif[1:3].upper() == "CG":
                type = "NCG"
            else:
                continue
        else:
            continue
        # Calculate GC content
        sequence = samfile.fetch(chrom, read.pos, read.pos + read.template_length).upper()
        gc_content = round((sequence.count('G') + sequence.count('C')) / len(sequence) *100)
        # Get ID
        id = str(gc_content) + '_' + str(read.template_length)
        # Add 'aberrant' and 'total' counts
        if id in outDict["ID"]:
            if type == 'CGN':
                outDict["CGN"][outDict["ID"].index(id)] += 1
            if type == 'NCG':
                outDict["NCG"][outDict["ID"].index(id)] += 1
        else:
            outDict["ID"].append(id)
            if type == 'CGN':
                outDict["CGN"].append(1)
                outDict["NCG"].append(0)
            if type == "NCG":
                outDict["NCG"].append(1)
                outDict["CGN"].append(0)
    # Disconnect .bam file from pysam library
    bamfile.close()
    # Return <outDict> as data.frame
    return pd.DataFrame.from_dict(outDict)

def fetchData(args, chromosome):
    # Define chromosome, start and end position
    chrom = chromosome.split('_')[0]
    start = int(chromosome.split('_')[1])
    end = int(chromosome.split('_')[2])
    # Load regions file
    if not args.type == "Whole_genome":
        regions = np.load(os.path.join(os.path.join(os.path.dirname(args.regions), args.type), f"{chrom}_{start}_{end}.npz"))['my_array']
    # Initiate output dictionary
    outDict = {"ID": [],
               "CGN": [],
               "NCG": []}
    # Load .sam file
    samfile = pysam.FastaFile(args.ref)
    # Conntect .bam file to pysam library
    bamfile = pysam.AlignmentFile(args.infile, "rb")
    # Extract relevant information from .bam file
    for read in bamfile.fetch(chrom, start, end):
        # Skip read that are not paired, properly paired, have low MAPQ value and are hard-clipped
        if not read.is_paired or not read.is_proper_pair or read.mapping_quality < 10 or 'H' in read.cigarstring or read.is_reverse:
            continue
        # Extract 5' end motif
        motif = read.query_alignment_sequence[:3]
        # Check if the end motif equals CGN or NCG
        if motif[:2].upper() == "CG":
            type = "CGN"
        elif motif[1:3].upper() == "CG":
            type = "NCG"
        else:
            continue
        if not args.type == "Whole_genome":
            # Check if 5' end falls within any region
            fallsWithin = (read.pos > regions[:, 0]) & (read.pos < regions[:, 1])
            if not fallsWithin.any():
                continue
        # Calculate GC content
        sequence = samfile.fetch(chrom, read.pos, read.pos + read.template_length).upper()
        gc_content = round((sequence.count('G') + sequence.count('C')) / len(sequence) *100)
        # Get ID
        id = str(gc_content) + '_' + str(read.template_length)
        # Add 'aberrant' and 'total' counts
        if id in outDict["ID"]:
            if type == 'CGN':
                outDict["CGN"][outDict["ID"].index(id)] += 1
            if type == 'NCG':
                outDict["NCG"][outDict["ID"].index(id)] += 1
        else:
            outDict["ID"].append(id)
            if type == 'CGN':
                outDict["CGN"].append(1)
                outDict["NCG"].append(0)
            if type == "NCG":
                outDict["NCG"].append(1)
                outDict["CGN"].append(0)
    # Disconnect .bam file from pysam library
    bamfile.close()
    # Return <outDict> as data.frame
    return pd.DataFrame.from_dict(outDict)


#-----------
# M A I N   S C R I P T
#-----------
def Main():
    args = ParsingArguments()
    # Transforms <contigs> into a list; splitting after each comma
    if args.type == "Tumor_specific":
        chr = args.contigs.split(',')
    else:
        chr = threading_parts(args=args, chunk=25000000)
    # Subset regions
    if not args.type == "Whole_genome" and not args.type == "Tumor_specific":
        subset_rpr(args=args, chromosomes=chr)
    # Create empty data.frame to populate
    output_df = pd.DataFrame()
    # Initiate threading
    pool = Pool(processes=args.threads)
    # Extract relevant information from .bam file
    if args.type == "Tumor_specific":
        fetch_partial = partial(fetchData_tumorSpecific, args)
        output_df = output_df.append(pool.map(fetch_partial, chr), ignore_index=True)
    else:  
        fetch_partial = partial(fetchData, args)
        output_df = output_df.append(pool.map(fetch_partial, chr), ignore_index=True)
    # Close pool
    pool.close()
    pool.join()
    # Sum 'aberrant' and 'total' fragments
    CGN = output_df.groupby('ID')['CGN'].sum().reset_index(name='CGN')
    NCG = output_df.groupby('ID')['NCG'].sum().reset_index(name='NCG')
    #non_aberrant = output_df.groupby('ID').size().reset_index(name='Non_aberrant')
    merged = pd.merge(CGN, NCG, on='ID', how='inner')
    # Add columns
    merged["GC"] = merged['ID'].str.split('_', expand=True)[0].astype(int)
    merged["Length"] = merged['ID'].str.split('_', expand=True)[1].astype(int)
    merged = merged.drop('ID', axis=1)
    # Re-order columns
    cols = merged.columns.tolist()
    cols = cols[2:] + [cols[0]] + [cols[1]]
    merged = merged[cols]
    # Re-order rows
    merged = merged.sort_values(by=['GC', 'Length'], ascending=[True, True])
    # Save <merged>
    merged.to_csv(args.outfile, sep='\t', index=False)

if __name__ == "__main__":
    Main()
