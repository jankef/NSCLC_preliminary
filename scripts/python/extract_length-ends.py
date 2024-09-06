#!/usr/bin/env python3
#
# Script Name: extract_length-ends.py
#
# Author: Florian Janke
# Last updated: 2024/09/06
#
# Description
#   (1) Extracts fragment end motives and fragment length of all fragments
#   (2) Summarizes count per GC content
#
#---------------------------------------------------------------------------------------------------

#-----------
# P A C K A G E S
#-----------
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
    parser.add_argument("-b", "--input",
                        dest="bamIn",
                        type=str,
                        required=True,
                        help="Path to the input bam file.")
    parser.add_argument("-n", "--nrBase",
                        dest="nrBase",
                        type=int,
                        default=3,
                        required=True,
                        help="Number of bases to fetch.")
    parser.add_argument("-m", "--mapq",
                        dest="mapq",
                        type=int,
                        required=True,
                        help="Applies MAPQ filter. Only to frequency table.")
    parser.add_argument("-o", "--output",
                        dest="output",
                        type=str,
                        required=True,
                        help="Path to output .tsv file.")
    parser.add_argument("-c", "--contigs",
                        dest="contigs",
                        type=str,
                        required=True,
                        help="List of contig names separated by a coma.")
    parser.add_argument("-t", "--threads",
                        dest="threads",
                        type=int,
                        required=True,
                        help="Number of threads.")
    parser.add_argument("-s", "--min_size",
                        dest="size_min",
                        default=0,
                        type=int,
                        required=True,
                        help="Minimum fragment length to consider a read.")
    parser.add_argument("-l", "--max_size",
                        dest="size_max",
                        default=250,
                        type=int,
                        required=True,
                        help="Maximum fragment length to consider a read.")
    parser.add_argument("-r", "--reference",
                        dest="ref",
                        type=str,
                        required=True,
                        help="Path to genome reference .fasta file.")
    return parser.parse_args()


#-----------
# F U N C T I O N S
#-----------
def getReadSeq(read, mate, nrBase):
    # Extract <nrBase> from both sides of the fragment
    if not read.is_reverse:
        P1_seq = read.query_alignment_sequence[:nrBase]
        P2_seq = mate.query_alignment_sequence[-nrBase:][::-1]  # Rev sequence.
    else:
        P2_seq = read.query_alignment_sequence[-nrBase:][::-1]  # Rev sequence.
        P1_seq = mate.query_alignment_sequence[:nrBase]
    # Return results
    return {"P1_seq": P1_seq,
            "P2_seq": P2_seq}

def CastDataTypes(data):
    # Cast object data into categorical,
    # and int32 for better memory usage and faster calculations.
    # print(data.info(memory_usage="deep"))
    data = data.reset_index(drop=True)
    for c in data.columns:
        if type(data[c][0]) == "object":
            data[c] = data[c].astype("category")
        elif type(data[c][0]) == str:
            data[c] = data[c].astype("category")
        elif type(data[c][0]) == "NA":
            data[c] = data[c].astype("category")
        elif type(data[c][0]) == int:
            data[c] = data[c].astype("int32")
        elif type(data[c][0]) == float:
            data[c] = data[c].astype("float32")
    return data

def fetchData(bamIn, nrBase, ref, chromosome):
    # Set the number of bases to look at; default = 3
    if nrBase:
        nrBase = nrBase
    else:
        nrBase = 3
    # Initiate output dictionary
    outDict = {"qname": [],
               "chr": [],
               "read_len": [],
               "strand": [],
               "start": [],
               "end": [],
               "left_seq": [],
               "right_seq": [],
               "mapq": [],
               "GC": []}
    # Load .sam file
    samfile = pysam.FastaFile(ref)
    # Conntect .bam file to pysam library
    bamfile = pysam.AlignmentFile(bamIn, "rb")
    readCache = {}
    # Extract relevant information from .bam file
    for read in bamfile.fetch(chromosome):
        # Filter hard-clipped reads as sequence ends might be clipped
        if 'H' in read.cigarstring:
            continue
        # Limit cache size 
        if len(readCache) > 10000:
            readCache = {}
        # Check if read length is greater 0
        if read.template_length == 0:
            continue
        # Place <read> into <readCache>. If the corresponding mate is found, sequence is extracted
        if read.query_name not in readCache:
            readCache[read.query_name] = read
            continue
        else:
            # Find mate in <readCache>
            mate = readCache[read.query_name]
            # Extract sequence
            SeqDict = getReadSeq(read, mate, nrBase)
            # Determine start and end position
            if not read.is_reverse:
                start_pos = read.reference_start
                end_pos = read.reference_start + abs(read.template_length)
            else:
                start_pos = mate.reference_start
                end_pos = mate.reference_start + abs(read.template_length)
            # Calculate GC content
            sequence = samfile.fetch(chromosome, start_pos, end_pos)
            sequence = sequence.upper()
            gc_content = round((sequence.count('G') + sequence.count('C')) / len(sequence) *100)
            # Fill <outDict> if 'P1' and 'P2' exist
            if ((SeqDict.get("P1_seq") != "") & (SeqDict.get("P2_seq") != "")):
                    outDict["qname"].append(read.query_name)
                    outDict["chr"].append(bamfile.get_reference_name(read.reference_id))
                    outDict["read_len"].append(abs(read.template_length))
                    outDict["strand"].append(0)
                    outDict["start"].append(start_pos)
                    outDict["end"].append(end_pos)
                    outDict["left_seq"].append(SeqDict.get("P1_seq"))
                    outDict["right_seq"].append(SeqDict.get("P2_seq"))
                    outDict["mapq"].append(read.mapping_quality)
                    outDict["GC"].append(gc_content)
        # Delete current query name from <readCache>
        del readCache[read.query_name]
    # Disconnect .bam file from pysam library
    bamfile.close()
    # Cast data types if <outDict> is not empty
    if len(pd.DataFrame.from_dict(outDict)) > 0:
        return CastDataTypes(pd.DataFrame.from_dict(outDict))
    else:
        return pd.DataFrame.from_dict(outDict)


#-----------
# M A I N   S C R I P T
#-----------
def Main():
    args = ParsingArguments()
    # Transforms <contigs> into a list; splitting after each comma
    ChrL = args.contigs.split(",")
    # Initiates data.frame
    OutDf = pd.DataFrame()
    # For .bam files >8Gb threading is turned off
    if (os.stat(args.bamIn).st_size / (1024*1024)) > 8000:
        args.threads = 1
    # Initiate threading
    pool = Pool(processes=args.threads)
    # Extract relevant information from .bam file
    fetchData_partial = partial(fetchData, args.bamIn, args.nrBase, args.ref)
    OutDf = OutDf.append(pool.map(fetchData_partial, ChrL), ignore_index=True)
    # Close pool
    pool.close()
    pool.join()
    # Remove duplicate rows
    OutDf.drop_duplicates(subset=["chr", "read_len", "start", "end"],
                        ignore_index=True,
                        inplace=True)
    # Remove rows containing 'Ns'
    OutDf = OutDf[~OutDf.left_seq.str.contains("N")]
    OutDf = OutDf[~OutDf.right_seq.str.contains("N")]
    # Filter by fragment length
    OutDf = OutDf[OutDf.read_len.isin(range(args.size_min, args.size_max))]
    # Filter by MAPQ
    OutDf = OutDf[OutDf.mapq >= args.mapq]
    # Calculate abundance per fragment length and GC content
    ends = ("left_seq", "right_seq")
    for end in ends:
        # Add 'id' column
        OutDf['id'] = OutDf['read_len'].astype(str) + '_' + OutDf['GC'].astype(str) + '_' + OutDf[end].astype(str)
        # Add 'Count' column
        OutDf['Count'] = 1
        # Group by 'id' and sum 'Count'
        tmp = OutDf.groupby('id')['Count'].sum()
        tmp = tmp.reset_index()
        tmp.columns = ['id', 'Count']
        # Add columns
        tmp["Length"] = tmp['id'].str.split('_', expand=True)[0].astype(int)
        tmp["GC"] = tmp['id'].str.split('_', expand=True)[1].astype(int)
        tmp["Motif"] = tmp['id'].str.split('_', expand=True)[2].astype(str)
        tmp["End"] = end
        tmp = tmp.drop('id', axis=1)
        # Re-order columns
        cols = tmp.columns.tolist()
        cols = cols[1:] + [cols[0]]
        tmp = tmp[cols]
        # Re-order rows
        tmp = tmp.sort_values(by=['Length', 'GC'], ascending=[True, True])
        if end == "left_seq":
            table = tmp
        else:
            table = table.append(tmp, ignore_index=True)
    # Create output directories
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    # Save fragment-ends as .tsv file
    table.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    Main()


