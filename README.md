# blast-align-tree
A pipeline to identify homologs and perform phylogenetic analysis with local BLAST databases
 
## Introduction
One common comparative genomic analysis is to find similar genes and compare them using phylogenetics methods. As an alternative to using online tools for such analyses, researchers may wish to download genomes of interest for local BLAST and downstream analyses. Homolog curation, tree construction, header parsing, and visualization alongside other quantitative data (e.g. expression data) can give quick insights into a gene family of interest.

![](pipeline.jpg)

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
bash blast-align-tree.sh AT1G19220.1 TAIR10cds.fa -n 30 30 30 -dbs TAIR10cds.fa Pvul218cds.fa Vung469cds.fa -hdr gene: polypeptide= locus= 
```

You can iterate on the tree by calling visualize-tree.R separately
```
Rscript visualize-tree.R -e AT1G19220.1 -b ARFs_v2
```

This should populate a subfolder with blast outputs, alignments, and tree visualizations

## Adding genomes
Add coding sequences files in fasta format to the genomes subdirectory, then compile a local BLAST database using 

```makeblastdb -in GenomeCDS.fa -parse_seqids -dbtype nucl```

Take note of appropriate header to specify when calling the bash script

