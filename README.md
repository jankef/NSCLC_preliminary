# Title to be added

## Background
to be added

## Authors
to be added

## Description
to bed added

## Usage
### Extract cfDNA fragment lengths and end motifs
This python script extracts cfDNA fragment length and fragment end nucleotide sequences simultaneously. The script returns a .tsv file containing the number of fragments of a specific length (bp), GC content (%), and fragment end motif. Fragment end motifs are extract from both the 5'- (left_seq) and 3'-end (right_seq).
```
python /path/to/script/extract_length-ends.py [options]
```
Parameters to set:
```
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


to be added

## Dependencies
### Software
- R version 4.4.0 (or newer)

### R packages
- ggplot2
- ggprism
- colorspace
- data.table
- VennDiagram
- stringr
- dplyr
- Hmisc
- survival
- GGally
- corrplot
- ggcorrplot
- cowplot
- showtext
- sva
- rlang


## License
This project is licensed under the APACHE 2.0 License - see the LICENSE.md file for details.
