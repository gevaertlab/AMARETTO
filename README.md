[//]: # (TODO: Bioconductor support?)
[//]: # (TODO: Some examples)

<p align="center">
  <a href="https://github.com/gevaertlab/AMARETTO">
    <img height="150" src="https://github.com/gevaertlab/AMARETTO/blob/master/files/logo.jpg">
  </a>
  <h1 align="center"></h1>
</p>


[![License](https://img.shields.io/badge/license-GPL-yellow.svg)](https://opensource.org/licenses/GPL-2.0)
[![Version](https://img.shields.io/badge/version-0.99.1-lightgrey.svg)]()

Integrating an increasing number of available multi-omics cancer data remains one of the main challenges to improve our understanding of cancer. One of the main challenges is using multi-omics data for identifying novel cancer driver genes. We have developed an algorithm, called AMARETTO, that integrates copy number, DNA methylation and gene expression data to identify a set of driver genes by analyzing cancer samples and connects them to clusters of co-expressed genes, which we define as modules. We applied AMARETTO in a pancancer setting to identify cancer driver genes and their modules on multiple cancer sites. AMARETTO captures modules enriched in angiogenesis, cell cycle and EMT, and modules that accurately predict survival and molecular subtypes. This allows AMARETTO to identify novel cancer driver genes directing canonical cancer pathways.

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Running AMARETTO](#running-amaretto)
- [References](#references)
- [License](#license)

## Introduction

This repository contains the code accompanying the manuscript ["Module analysis captures pancancer genetically and epigenetically deregulated cancer driver genes for smoking and antiviral response"](https://www.sciencedirect.com/science/article/pii/S2352396417304723).
We have developed an algorithm called AMARETTO to identify pancancer driver genes. AMARETTO integrates pancancer DNA copy number, DNA methylation, and gene expression data into modules to identify pancancer driver genes. 

The algorithm: 
* Step 1 identifies candidate cancer driver genes with tumor-specific DNA copy number or DNA methylation alterations compared to normal tissue. We used GISTIC to identify significantly and recurrently deleted or amplified regions in the genome. Similarly, we used MethylMix to identify recurrently hyper- or hypomethylated genes. 
* Step 2 identifies cancer driver genes by modeling the relationship between (epi)genomic and transcriptomic data on an individual gene basis.
* Step 3 uses the cancer driver genes identified from Step 2 and takes a global approach by dissecting global gene expression data into modules of co-expressed genes. Each module also has an associated regulatory program that connects the cancer driver genes from Step 1 with their downstream targets. This regulatory program is modeled using linear regression with elastic net regularization.
 
AMARETTO supports downloading and processing TCGA data from [Firehose](https://gdac.broadinstitute.org/).

## Installation

Install from the GitHub repository using devtools:

    install.packages("devtools")
    library(devtools)
    devtools::install_github("gevaertlab/AMARETTO")

## Running AMARETTO

* The RUNNING_AMARETTO_EXAMPLE folder contains an example R script for a typical AMARETTO analysis. Please try!
* [*Introduction to AMARETTO*](https://nbviewer.jupyter.org/github/gevaertlab/AMARETTO/blob/master/vignettes/amaretto.pdf) provides a comprehensive example of the AMARETTO workflow with detailed  explanations of each function.</br> 
* Detailed information on `AMARETTO` package functions can be obtained in the help files. For example, to view the help file for the function `AMARETTO` in a R session, use `?AMARETTO`.

## References

1.	Champion, M. et al. Module Analysis Captures Pancancer Genetically and Epigenetically Deregulated Cancer Driver Genes for Smoking and Antiviral Response. EBioMedicine 27, 156–166 (2018).
2.	Gevaert, O., Villalobos, V., Sikic, B. I. & Plevritis, S. K. Identification of ovarian cancer driver genes by using module network integration of multi-omics data. Interface Focus 3, 20130013–20130013 (2013).
3.	Gevaert, O. MethylMix: an R package for identifying DNA methylation-driven genes. Bioinformatics 31, 1839–1841 (2015).


## License

Copyright (c) 2018 Jayendra Shinde

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
