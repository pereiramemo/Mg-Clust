# MG-Clust
## Clustering of ORF sequences in metagenomic data.
This repository contains the pipeline mg_clust. This tool is a command line application programmed in BASH and AWK, dedicated to the computation of Operational Protein Units (OPUs) in metagenomic data, based on Open Reading Frame (ORF) amino acid sequences. It takes as an input a sampling set of preprocessed (unassembled) metagenomic samples and outputs the OPUs abundance table.

Dependencies:  
[MEGAHIT](https://github.com/voutcn/megahit)  
[BWA](http://bio-bwa.sourceforge.net)  
[SAMTools](https://github.com/samtools/)  
BEDTools(https://bedtools.readthedocs.io/en/latest/)  
Picard(https://broadinstitute.github.io/picard/)  
MMSeqs2(https://github.com/soedinglab/MMseqs2)  
[GNU Parallel](https://www.gnu.org/software/parallel/)  

![MG-Clust workflow](./figures/MG-Clust.png)


