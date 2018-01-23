# bioinf_sparc
An implementation of the [Sparc](https://peerj.com/articles/2016.pdf) algorithm.
Given a backbone, reads and a mapping file describing where the reads align with the backbone, creates a sparse k-mer graph with weighted edges. Consensus is found by traversing the graph and finding the heaviest path from the beginning to the end of the backbone.

Saves output to "output.fasta".

Generate mapping file by using [minimap2](https://github.com/lh3/minimap2), be sure to add the -c tag to generate CIGAR data:
minimap2 -c path_to_backbone.fasta path_to_reads.fastq > path_to_mappings.paf

To compile:
g++ -std=c++11 -o Sparc Sparc.cpp

To run:
Sparc -g x -k y -b path_to_backbone.fasta -r path_to_reads.fastq  -m path_to_mappings.paf 

