# deNOPA

### Introduction

*A*s the basal bricks, the dynamics and arrangement of nucleosomes orchestrate the higher architecture of chromatin in a fundamental way, thereby affecting almost all nuclear biology processes. Thanks to its rather simple protocol, ATAC-seq has been rapidly adopted as a major tool for chromatin-accessible profiling at both bulk and single-cell level; however, to picture the arrangement of nucleosomes *per se* remains a challenge with ATAC-seq. We introduce a novel ATAC-seq analysis toolkit, named **deNOPA**, to predict nucleosome positions. Assessments showed that deNOPA not only outperformed state-of-the-art tools, but it is the only tool able to predict nucleosome position precisely with ultrasparse ATAC-seq data. The remarkable performance of deNOPA was fueled by the reads from short fragments, which compose nearly half of sequenced reads and are normally discarded from an ATAC-seq library. However, we found that the short fragment reads enrich information on nucleosome positions and that the linker regions were predicted by reads from both short and long fragments using Gaussian smoothing. 

### Install

##### Pre-requirements

The deNOPA package was initially developed using python 2.7. The support of python 3 has also been added. The package was tested under the default environment of Anaconda-5.3.1, both python 2.7 and python 3.7 version. Besides a python environment, the following dependencies were also needed. 

* numpy
* scipy
* h5py
* pysam
* sklearn
* statsmodels

Please make sure they were properly installed ahead of the deNOPA package itself. 

##### Install from source code

Use the following commands to get deNOPA installed. 

```
git clone https://gitee.com/bxxu/denopa.git

cd denopa

python setup.py install
```

##### Install from pre-built files. 

Download the compatible wheel file from the dist directory in this repo according to your version of python. Then get it installed using the following command:

```
pip install deNOPA-x.y.z-pyX-none-any.whl
```

### Usage

```
usage: denopa [-h] -i INPUT [-o OUTPUT] [-b BUFFERSIZE] [-s CHROMSKIP] [-n NAME] [-m MAXLEN] [--proc PROC] [-p PARER] [-q QNFR] [-r]

Decoding the nucleosome positions with ATAC-seq data at single cell level

optional arguments:

-h, --help            show this help message and exit

-i INPUT, --input INPUT

​                         The input bam files. The files should be sorted. This argument could be given multiple times for multiple input files.

-o OUTPUT, --output OUTPUT

​                         The directory where the output files will go to. It will be created if not exists (default .).

-b BUFFERSIZE, --bufferSize BUFFERSIZE

​                         Number of reads buffered in reading the bam file (default 1000000).

-s CHROMSKIP, --chromSkip CHROMSKIP

​                         Names of chromosomes skiped from the processing. Multiple values should be sepaated by ',' (default chrY,chrM).

-n NAME, --name NAME  The name of the project (default deNOPA).

-m MAXLEN, --maxLen MAXLEN

​                         The maximun fragment length in the input files (default 2000).

--proc PROC

​                         Number of processors used in the analysis (default 1).

-p PARER, --pARER PARER

​                         The p-value threshold used in determining the ATAC-seq reads enriched regions (ARERs, default 0.1).

-q QNFR, --qNFR QNFR  

​                         The q-value threshold used in determining the nucleosome free regions (NFRs, default 0.1).

-r, --removeIntermediateFiles

​                         The intermediate files will be removed if this flag is set.
```

### Output files

##### Results

* {NAME}_ARERs.txt: The bet like file of the detected ARERs. The chromosome name, start position, end position, local maximum point, maximum signal and p-value of each ARER were recorded. 
* {NAME}_NFR.txt: A standard broadPeak formatted bed file for detected NFRs. 
* {NAME}_nucleosomes.txt: A bed like file for all detected nucleosomes. For each nucleosome, the chromosome name, start position, end position, left inflection point, center position, right inflection point, number of sequencing fragments intersected, number of sequencing fragments covered it, number of Tn5 cutting sites in its inner, p-value indicating whether the number of sequencing fragments covering the nucleosome was large enough, p-value indicating whether the number of Tn5 cutting events in its inner was small enough, the combination of the about two p-values, whether  the nucleosome was dynamic were listed respectively. 

##### Intermediated files

* {NAME}_candidates.pkl: All the detected nucleosome candidates before the DBSCAN outlier detection was applied. 
* {NAME}_frag_len.pkl: The fragment length distribution of the ATAC-seq library. 
* {NAME}_pileup_signal.hdf: The raw coverage and cutting sites signal profiles. 
* {NAME}_smooth.hdf: The smoothed signal profiles and their derivations. 

### Citation

Please add the following citation if you use deNOPA in your study: 

> Xu B, Li X, Gao X, et al. DeNOPA: decoding nucleosome positions sensitively with sparse ATAC-seq data[J]. Briefings in Bioinformatics, 2022, 23(1): bbab469.

### Contacts

You can send questions, discussions, bug reports and other useful information to Zhihua Zhang (zangzhihua@big.ac.cn) or Bingxiang Xu (xubingxiang@sus.edu.cn).























