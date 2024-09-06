# Title to be added

## Background
to be added

## Authors
to be added

## Description
The code provided in this repository allows the extraction of the cfDNA fragmentation features evaluated in 'PAPER TITLE'. These features include: (i) cfDNA fragment length, (ii) cfDNA fragment end trinucleotides, as well as (iii) aberrant cfDNA fragment end positions. In addition, motif diversity scores (MDS), Alu CGN/NCG ratios, and information-weighted fraction of aberrant fragments (iwFAF) scores is provided. A reduced dataset to test the provided scripts can be found in /example_data/. The .bed file containing the location of recurrently protected regions was downloaded from Budhraja 2023 (doi:10.1126/scitranslmed.abm6863) and lifted to hg38. Alu element regions were downloaded from the UCSC Table Browser. Parts of the provided code were apadated from Moldovan 2024 (https://doi.org/10.1016/j.xcrm.2023.101349). The calculation of P126-135 and D126-135 are detailed in the Supplemenatry Methods and not described here. The calculation of ctCPA scores was perviously described (Janke 2024,  https://doi.org/10.1002/ijc.35152) and can be found here: https://github.com/jankef/ctDNA-predicts-outcome-in-head-and-neck-cancer.


## Manual
### cfDNA fragmentation feature extraction
#### Fragment length and end motifs
This python script extracts cfDNA fragment length and fragment end nucleotide sequences simultaneously. The script returns a .tsv file containing the number of fragments of a specific length (bp), GC content (%), and fragment end motif. Fragment end motifs are extract from both the 5'- (left_seq) and 3'-end (right_seq).

```
python /path/to/script/extract_length-ends.py [options]
```

##### Parameters to set:
```
Options:
--input     [Path to the input .bam file]
--nrBase    [Number of bases to extract from each read; Default = 3]
--mapq      [Minimum MAPQ value of a read to be considered for the extraction of fragment length and end motif]
--output    [Path to .tsv output file]
--contigs   [Comma separated list of chromosomes to be considered. E.g., 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22']
--threads   [Number of threads to be used]
--min_size  [Minimum fragment length to consider a read for the extraction of fragment length and end motif; Default = 0]
--max_size  [Maximum fragment length to consider a read for the extraction of fragment length and end motif; Default = 250]
--reference [Path to reference genome .fasta file]
```

##### Example usage
```
python /path/to/script/extract_length-ends.py \
  --input /path/to/example.bam \
  --nrBase 3 \
  --mapq 30 \
  --output /path/to/output.tsv \
  --contigs 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22' \
  --threads 10 \
  --min_size 0 \
  --max_size 250 \
  --reference /path/to/reference.fasta
```

#### Aberrant cfDNA fragment end positions
Counts the number of aberrant and non-aberrant cfDNA fragments per GC content (%) and fragment length (bp). Returns a .tsv file with that can be used to calculate iwFAF scores.

```
python /path/to/script/aberrant_positions.py [options]
```

##### Parameters to set:
```
Options:
--infile     [Path to the input .bam file]
--outfile    [Path to .tsv output file]
--rpr        [Path to the .bed file containing the positions of recurrently protected regions.]
--contigs    [Comma separated list of chromosomes to be considered. E.g., 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22']
--reference  [Path to reference genome .fasta file]
--threads    [Number of threads to be used]
```

##### Example usage
```
python /path/to/script/aberrant_positions.py \
    --infile /path/to/example.bam \
    --outfile /path/to/output.tsv \
    --contigs 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22' \
    --threads 10 \
    --reference /path/to/reference.fasta \
    --rpr /path/to/rpr_hg38.bed
```

#### CGN/NCG ratios at Alu element regions
Counts the number of cfDNA fragments that carry a CGN or NCG end motif, separately. The number of end motif occurrences are returned per fragment length (bp) and GC content (%). Returns a .tsv file.

```
python /path/to/script/alu_element-ratio.py [options]
```

##### Parameters to set:
```
Options:
--infile     [Path to the input .bam file]
--outfile    [Path to .tsv output file]
--regions    [Path to the .bed file containing the positions of alu elements]
--contigs    [Comma separated list of chromosomes to be considered. E.g., 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22']
--reference  [Path to reference genome .fasta file]
--threads    [Number of threads to be used]
--type       [Specifies what kind of .bed file was provided. Accepts 'Alu', 'Whole_genome', and 'CpG_islands']
```

##### Example usage
```
python /path/to/script/alu_element-ratio.py \
    --infile /path/to/example.bam \
    --outfile /path/to/output.tsv \
    --contigs 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22' \
    --threads 10 \
    --reference /path/to/reference.fasta \
    --regions /path/to/alu_element_hg38.bed \
    --type 'Alu'
```

### Biomarker calculation
#### Motif diversity score
This R script calculates iwFAF scores from the output of 'extract_length-ends.py'. All input .tsv files need to be in the same directory.

```
Rscript /path/to/script/MDS_calc.R [options]
```

##### Parameters to set:
```
Options:
--input     [Path to the directory containing the 'extract_length-ends.py' output .tsv files.]
--output    [Path to .tsv output file]
--batch     [Path to a .tsv file specifying the batch of the individual samples in 'input'. This file has two columns: Sample_ID and Batch_ID. Batch_ID is an integer that specifies which sample belongs to batch 1, 2, 3, etc.]
```


#### iwFAF score
This R script calculates iwFAF scores from the output of 'aberrant_positions.py'. All input .tsv files need to be in the same directory. File names starting with 'CTRL' are considered as controls and used for calculating the probability of a fragment being normal/aberrant. When iwFAF scores are calculated for controls. All expect for the currently assessed controls are used for probability calculation.

```
Rscript /path/to/script/iwFAF-score_calc.R [options]
```

##### Parameters to set:
```
Options:
--input     [Path to the directory containing the 'aberrant_positions.py' output .tsv files.]
--output    [Path to .tsv output file]
```


to be added

## Dependencies
### Software
- R version 4.4.0 (or newer)
- Python version 3.7.3 (or newer)

### R packages
- data.table
- matrixStats
- stringr
- optparse

### Python modules
- os
- argparse
- multiprocessing
- functools
- Pandas
- pysam
- Bio.SeqUtils
- FrEIA_tools
- numpy

## License
This project is licensed under the APACHE 2.0 License - see the LICENSE.md file for details.
