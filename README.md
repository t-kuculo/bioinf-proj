# bioinf_sparc
An implementation of the [Sparc](https://peerj.com/articles/2016.pdf) algorithm.
Given a backbone, reads and a mapping file describing where the reads align with the backbone, creates a sparse k-mer graph with weighted edges. Consensus is found by traversing the graph and finding the heaviest path from the beginning to the end of the backbone.

Saves output to "output.fasta".

To compile:
g++ -std=c++11 -o Sparc Sparc.cpp

To run:
Sparc -g x -k y -b path_to_backbone.fasta -r path_to_reads.fastq  -m path_to_mappings.paf 
