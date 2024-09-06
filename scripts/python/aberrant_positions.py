#!/usr/bin/env python3
#
# Script Name: aberrant_positions.py
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
    parser.add_argument("-x", "--rpr",
                        dest="rpr",
                        type=str,
                        required=True,
                        help="Path to recurrently protected regions .bed file.")
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

def get_status(data, start, end):
    data = data[end - data[:, 0] > 0]
    if len(data) == 0:
            return 'exclude'
    data = data[(data[:, 1] - start) > 0]
    if len(data) == 0:
            return 'exclude'
    array = (start > data[:, 0]) & (start < data[:, 1])
    if array.any():
        return 'aberrant'
    array2 = (end > data[:, 0]) & (end < data[:, 1])
    if array2.any():
        return 'aberrant'
    else:
        return 'non-aberrant'

def fetchData(args, chromosome):
    # Define chromosome, start and end position
    chrom = chromosome.split('_')[0]
    start = int(chromosome.split('_')[1])
    end = int(chromosome.split('_')[2])
    # Load RPR file
    rpr = np.load(os.path.join(os.path.join(os.path.dirname(args.rpr), 'chunks'), f"{chrom}_{start}_{end}.npz"))['my_array']
    # Initiate output dictionary
    outDict = {"ID": [],
               "aberrant": [],
               "non_aberrant": []}
    # Load .sam file
    samfile = pysam.FastaFile(args.ref)
    # Conntect .bam file to pysam library
    bamfile = pysam.AlignmentFile(args.infile, "rb")
    # Extract relevant information from .bam file
    for read in bamfile.fetch(chrom, start, end):
        # Skip read that are not paired, properly paired, have low MAPQ value and are hard-clipped
        if not read.is_paired or not read.is_proper_pair or read.mapping_quality < 10 or 'H' in read.cigarstring or read.is_reverse:
            continue
        # Get 'aberrant', 'non-aberrant' and 'exclude' status
        status = get_status(data=rpr, start=read.pos, end=read.pos + read.template_length)
        if status == 'exclude':
            continue
        # Calculate GC content
        sequence = samfile.fetch(chrom, read.pos, read.pos + read.template_length).upper()
        gc_content = round((sequence.count('G') + sequence.count('C')) / len(sequence) *100)
        # Get ID
        id = str(gc_content) + '_' + str(read.template_length)
        # Add 'aberrant' and 'total' counts
        if id in outDict["ID"]:
            if status == 'aberrant':
                outDict["aberrant"][outDict["ID"].index(id)] += 1
            if status == 'non-aberrant':
                outDict["non_aberrant"][outDict["ID"].index(id)] += 1
        else:
            outDict["ID"].append(id)
            if status == 'aberrant':
                outDict["aberrant"].append(1)
                outDict["non_aberrant"].append(0)
            if status == "non-aberrant":
                outDict["non_aberrant"].append(1)
                outDict["aberrant"].append(0)
    # Disconnect .bam file from pysam library
    bamfile.close()
    # Cast data types if <outDict> is not empty
    if len(pd.DataFrame.from_dict(outDict)) > 0:
        return CastDataTypes(pd.DataFrame.from_dict(outDict))
    else:
        return pd.DataFrame.from_dict(outDict)

def subset_rpr(args, chromosomes):
    # Create output folder
    chunk_dir = os.path.join(os.path.dirname(args.rpr), 'chunks')
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
    # Load RPR file
    rpr = pd.read_csv(args.rpr, sep="\t", header=None, usecols=[0, 1, 2])
    rpr.columns = ["chr", "start", "end"]
    # Calculate RPR width
    rpr['width'] = rpr['end'] - rpr['start']
    for chunk in chromosomes:
        chunk = chunk.split('_')
        # Skip if file already existis
        if os.path.exists(os.path.join(chunk_dir, f"{chunk[0]}_{chunk[1]}_{chunk[2]}.npz")):
            continue
        # Subset
        tmp = rpr[rpr['chr'] == chunk[0]]
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


#-----------
# M A I N   S C R I P T
#-----------
def Main():
    args = ParsingArguments()
    # Transforms <contigs> into a list; splitting after each comma
    chr = threading_parts(args=args, chunk=25000000)
    # Subset RPRs
    subset_rpr(args=args, chromosomes=chr)
    # Create empty data.frame to populate
    output_df = pd.DataFrame()
    # Initiate threading
    pool = Pool(processes=args.threads)
    # Extract relevant information from .bam file
    fetch_partial = partial(fetchData, args)
    output_df = output_df.append(pool.map(fetch_partial, chr), ignore_index=True)
    # Close pool
    pool.close()
    pool.join()
    # Sum 'aberrant' and 'total' fragments
    aberrant = output_df.groupby('ID')['aberrant'].sum().reset_index(name='Aberrant')
    non_aberrant = output_df.groupby('ID')['non_aberrant'].sum().reset_index(name='Non_aberrant')
    #non_aberrant = output_df.groupby('ID').size().reset_index(name='Non_aberrant')
    merged = pd.merge(non_aberrant, aberrant, on='ID', how='inner')
    # Add columns
    merged["GC"] = merged['ID'].str.split('_', expand=True)[0].astype(int)
    merged["Length"] = merged['ID'].str.split('_', expand=True)[1].astype(int)
    merged = merged.drop('ID', axis=1)
    # Re-order columns
    cols = merged.columns.tolist()
    cols = cols[2:] + [cols[1]] + [cols[0]]
    merged = merged[cols]
    # Re-order rows
    merged = merged.sort_values(by=['GC', 'Length'], ascending=[True, True])
    # Save <merged>
    merged.to_csv(args.outfile, sep='\t', index=False)

if __name__ == "__main__":
    Main()