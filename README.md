# demultiplexer for nanopore

Chances are that if you have been doing some metabarcoding, you have done it with Illumina platforms. But now you want to try nanopore. Oxford nanopore offers its propietary barcoding indices, which allows you to demultiplex up to 96 samples per run. Which is probably all you need in most cases. BUT! What if you want to pool 1000s of samples per run of a flowcell?

Here I have modified a Illumina pipeline for custom indices (github.com/ramongallego/demultiplexer_for_dada2), which I did following the backbone of **banzai** (github.com/jimmyodonnell/banzai).

The idea is that you have a `metadata` file with the information we need: sample name, index at the 5' end (which I named i5 because it is hard to move from Illumina), index at the 3' end and so on. Then you need a parameters file, in which you point at the right files, and establish some settings: minimum and maximum fragment length, number of cores (checked in Linux, let me know if it doesn´t work in your OS) and so on.

The data I had was generated using 9.4 flowcells, and although I used `sup`basecalling, there are still some errors that can lead to low efficiency of demultiplexing. So I used *anchors*: sequences immediately preceding the indices so we can tell `cutadapt` where exactly they are. Once I get my hands on some kit14 and last generation flowcells, this probably won´t be neccessary.

## Dependencies

This script has been tested in Pop OS 21. It should work in Mac as well, but I think different `sed` versions will return different results.

The script relies mostly on cutadapt (<https://cutadapt.readthedocs.io/en/stable/index.html>).

## Clustering

There are not that many algorithms out there for clustering Nanopore-generated amplicons. `decona` is one, and to run it you need a Linux machine with `conda`, and activate the conda environment. So my approach is to run the demultiplexer first, activate the environment and then run `bash decona.sh <path/to/output/folder>`pointing to the outputfolder you just generated. This will add all decona's output folders, and also two new output files:

  * `ASV_table.csv` is the occurrence table, in long format. It has the columns: Sample, Hash, nReads
  * `Hash_key.csv` has the DNA sequences of each Hash
  
