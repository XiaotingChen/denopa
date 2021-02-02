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

> 

