# ZikaPhyloSurveillanceWGS
#### Python code for analyses performed in the Zika virus whole genome sequencing manuscript published in Nature Scientific Reports, "Phylogenetic surveillance of travel-related Zika virus infections through whole-genome sequencing methods"

## Usage

### Patristic distance analysis

The following script accepts phylogenetic tree in newick format as input.
The tips should be labelled in the following format in order to identify which country/group each belongs to:

'>1_Country'

python PatDistSpectrum.v1.py -i infile -o1 outfile1 -o2 outfile2

Outfile1 being the total patristic distances in the phylogenetic tree.

Outfile2 being the intra-country comparisons using patristic distance percentiles. 

In this study only the min (0th percentile) and the max (100th percentile) patristic distance  were utilized.

### Unique variant analysis

Nucleotide counting script
This accepts an aligned fasta file as input, again requiring the ids to be labelled as illustrated above.
It outputs a file to be used in the following script, varianceDetection.v1.py
The specified output file is a csv file with all counts of the different nucleotides at each position.
It can process sequences where there are only 1 per country, or, multiple seqs per country (specified with the -t flag).

python nucFreqAlign.v1.py -h
usage: nucFreqAlign.v1.py [-h] -i INFILE -t TYPE -r REDUNDANT -o OUTFILE

optional arguments:
  -h, --help            show this help message and exit
  
  -i INFILE, --infile INFILE
                        Specify path to input alignment file in fasta format
                        
  -t TYPE, --type TYPE  Specify file type: multiple sequences per group/country? '-t m', or single sequence per group? '-t s'
                        
  -r REDUNDANT, --redundant REDUNDANT
                        Specify whether redundant bases are to be incorporated                        
                        in analysis (T) or ignored (F) or '
                        
  -o OUTFILE, --outfile OUTFILE
                          Specify output filename

This accepts input from the nucFreqAlign.v1.py file

It outputs two files:

1) A NucleotideVariance file which is similar to the input file except that it flags positions that have inter-country variance, ie, different between countries.

2) A CountryVariance file which  lists each variant as well as how many countries share it.

If multiple sequences are present per country, the frequencies for each country are also listed.

python varianceDetection.v1.py -h

usage: varianceDetection.v1.py [-h] -i INFILE -t TYPE

optional arguments:
  -h, --help            show this help message and exit
  
  -i INFILE, --infile INFILE
Specify path to input csv file
                        
  -t TYPE, --type TYPE  Specify file type: multiple sequences per
group/country? '-t m', or single sequence per group? '-t s'

### Shannon diversities from deep sequencing data

The following script outputs a csv file with counts for each nucleotide as well as the shannon diversity at each position.
The input file is an mpileup file generated in samtools.

python mpileup2snp.v1.py -h

usage: mpileup2snp.v1.py [-h] -i MPILEUP -d DEPTH -b STRANDBIAS -r REFERENCELENGTH

optional arguments:
  -h, --help            show this help message and exit
  
  -i MPILEUP, --mpileup MPILEUP Specify path to input mpileup file
                        
  -d DEPTH, --depth DEPTH
  
                        Input a depth threshold for variant calling (ex: 0.1 where 0.1 = 10 percent depth). If not important, enter '0.0'
                        
  -b STRANDBIAS, --strandbias STRANDBIAS
  
                        Input a threshold if strand bias is required (ex 100, where at least a 100 for and rev reads are required). If not important, enter '0'
                        
  -r REFERENCELENGTH, --referenceLength REFERENCELENGTH
                          Input a reference length, if multiple segments are present, use the longest reference length
