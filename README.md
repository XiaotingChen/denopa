# deNOPA

### Introduction

*A*s the basal bricks, the dynamics and arrangement of nucleosomes orchestrate the higher architecture of chromatin in a fundamental way, thereby affecting almost all nuclear biology processes. Thanks to its rather simple protocol, ATAC-seq has been rapidly adopted as a major tool for chromatin-accessible profiling at both bulk and single-cell level; however, to picture the arrangement of nucleosomes *per se* remains a challenge with ATAC-seq. We introduce a novel ATAC-seq analysis toolkit, named **deNOPA**, to predict nucleosome positions. Assessments showed that deNOPA not only outperformed state-of-the-art tools, but it is the only tool able to predict nucleosome position precisely with ultrasparse ATAC-seq data. The remarkable performance of deNOPA was fueled by the reads from short fragments, which compose nearly half of sequenced reads and are normally discarded from an ATAC-seq library. However, we found that the short fragment reads enrich information on nucleosome positions and that the linker regions were predicted by reads from both short and long fragments using Gaussian smoothing. 

### Install

##### Pre-requirements

The current version of deNOPA was developed with python 2.7. Python 3 is not supported by now. Besides a python environment, you need also get the following pre-requirements installed:

* numpy
* scipy
* h5py
* pysam
* sklearn
* statsmodels

##### Installation

Use the following commands to get deNOPA installed. 

> git clone https://gitee.com/bxxu/denopa.git
>
> cd denopa
>
> python setup.py install

### Usage

> usage: denopa [-h] -i INPUT [-o OUTPUT] [-b BUFFERSIZE] [-s CHROMSKIP]
>               [-n NAME] [-m MAXLEN] [--proc PROC] [-p PARER] [-q QNFR] [-r]
>
> Decoding the nucleosome positions with ATAC-seq data at single cell level
>
> optional arguments:
>   -h, --help            show this help message and exit
>   -i INPUT, --input INPUT
>                         The input bam files. The files should be sorted. This
>                         argument could be given multiple times for multiple
>                         input files.
>   -o OUTPUT, --output OUTPUT
>                         The directory where the output files will go to. It
>                         will be created if not exists (default .).
>   -b BUFFERSIZE, --bufferSize BUFFERSIZE
>                         Number of reads buffered in reading the bam file
>                         (default 1000000).
>   -s CHROMSKIP, --chromSkip CHROMSKIP
>                         Names of chromosomes skiped from the processing.
>                         Multiple values should be sepaated by ',' (default
>                         chrY,chrM).
>   -n NAME, --name NAME  The name of the project (default deNOPA).
>   -m MAXLEN, --maxLen MAXLEN
>                         The maximun fragment length in the input files
>                         (default 2000).
>   --proc PROC           Number of processors used in the analysis (default 1).
>   -p PARER, --pARER PARER
>                         The p-value threshold used in determining the ATAC-seq
>                         reads enriched regions (ARERs, default 0.1)
>   -q QNFR, --qNFR QNFR  The q-value threshold used in determining the
>                         nucleosome free regions (NFRs, default 0.1).
>   -r, --removeIntermediateFiles
>                         The intermediate files will be removed if this flag is
>                         set.

























