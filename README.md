# blast-align-tree
 A pipeline to visualize homologs in locally stored genomes
 
## Introduction
A simple comparative genomic analysis is to find similar genes and compare them using phylogenetics methods. Given the pace of change in plant genomic resources, it is helpful to compile a local set of genomes with associated BLAST databases. Simple BLAST searches followed by tree building, header parsing, and tree visualization can give a rapid view of a gene family of interest.

(pipeline picture to come)

## Installation
The following should be in your path

clustalo (ClustalOmega)

FastTree http://www.microbesonline.org/fasttree/

prank http://wasabiapp.org/software/prank/

trimAl http://trimal.cgenomics.org/trimal

Clone the repository which includes two genomes, Arabidopsis (TAIR10) and common bean 2.1

An example output can be found in subfolder AT4G33430.1
A markdown workbook of the tree visualization for AT4G33430.1 can be found in /examples

Try the pipline yourself by pulling 30 homologs of Arabidopsis ARF19 
```
sh blast-align-tree.sh AT1G19220.1 TAIR10cds.fa -n 30 30 30 -dbs TAIR10cds.fa Pvul218cds.fa Vung469cds.fa -hdr gene: polypeptide= locus= 
```

This should populate a subfolder with blast outputs, alignments, and tree visualizations

## Adding genomes
Add coding sequences files in fasta format to the genomes subdirectory, then compile a local BLAST database using "makeblastdb -in GenomeCDS.fa -parse_seqids -dbtype nucl". Take note of appropriate header to specify when calling the bash script

