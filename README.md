[![R build status](https://github.com/xihaoli/STAARpipelineSummary/workflows/R-CMD-check/badge.svg)](https://github.com/xihaoli/STAARpipelineSummary/actions)
[![Build Status](https://travis-ci.com/xihaoli/STAARpipelineSummary.svg?branch=main)](https://app.travis-ci.com/github/xihaoli/STAARpipelineSummary)
[![Build status](https://ci.appveyor.com/api/projects/status/glc569adbr7ar11j/branch/main?svg=true)](https://ci.appveyor.com/project/xihaoli/staarpipelinesummary/branch/main)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# STAARpipelineSummary
This is an R package for summarizing and visualizing association analysis results of whole-genome/whole-exome sequencing (WGS/WES) studies generated by STAARpipeline.
## Description
STAARpipelineSummary is an R package for summarizing and visualizing association analysis results generated by <a href="https://github.com/xihaoli/STAARpipeline">STAARpipeline</a>.
## Prerequisites
<a href="https://www.r-project.org">R</a> (recommended version >= 3.5.1)

For optimal computational performance, it is recommended to use an R version configured with the Intel Math Kernel Library (or other fast BLAS/LAPACK libraries). See the <a href="https://software.intel.com/en-us/articles/using-intel-mkl-with-r">instructions</a> on building R with Intel MKL.
## Dependencies
STAARpipelineSummary imports R packages <a href="https://cran.r-project.org/web/packages/Rcpp/index.html">Rcpp</a>, <a href="https://github.com/xihaoli/STAAR">STAAR</a>, <a href="https://github.com/xihaoli/STAARpipeline">STAARpipeline</a>, <a href="https://github.com/zilinli1988/SCANG">SCANG</a>, <a href="https://cran.r-project.org/web/packages/dplyr/index.html">dplyr</a>, <a href="https://bioconductor.org/packages/release/bioc/html/SeqArray.html">SeqArray</a>, <a href="https://bioconductor.org/packages/release/bioc/html/SeqVarTools.html">SeqVarTools</a>, <a href="https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html">GenomicFeatures</a>, <a href="https://bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg38.knownGene.html">TxDb.Hsapiens.UCSC.hg38.knownGene</a>, <a href="https://cran.r-project.org/web/packages/GMMAT/index.html">GMMAT</a>, <a href="https://bioconductor.org/packages/release/bioc/html/GENESIS.html">GENESIS</a>, <a href="https://cran.r-project.org/web/packages/Matrix/index.html">Matrix</a>, <a href="https://cran.r-project.org/web/packages/lattice/index.html">lattice</a>. These dependencies should be installed before installing STAARpipelineSummary.
## Installation
```
library(devtools)
devtools::install_github("xihaoli/STAARpipelineSummary",ref="main")
```
## Docker Image
A [docker image for STAARpipelineSummary](https://hub.docker.com/repository/docker/zilinli/staarpipeline), including R (version 3.6.1) built with Intel MKL and all STAAR-related packages (STAAR, SCANG, STAARpipeline, STAARpipelineSummary) pre-installed, is located in the Docker Hub. The docker image can be pulled using
```
docker pull zilinli/staarpipeline:0.9.6
```
## Usage
Please see the <a href="docs/STAARpipelineSummary_manual.pdf">**STAARpipelineSummary** user manual</a> for detailed usage of STAARpipelineSummary package. Please see the <a href="https://github.com/xihaoli/STAARpipeline-Tutorial">**STAARpipeline** tutorial</a> for a detailed example of summarizing and visualizing association analysis results using STAARpipelineSummary.
## Data Availability
The whole-genome individual functional annotation data assembled from a variety of sources and the computed annotation principal components are available at the [Functional Annotation of Variant - Online Resource (FAVOR)](https://favor.genohub.org) site.
## Version
The current version is 0.9.6 (April 16, 2022).
## Citation
If you use **STAARpipeline** and **STAARpipelineSummary** for your work, please cite:

Zilin Li*, Xihao Li*, Hufeng Zhou, Sheila M. Gaynor, Margaret S. Selvaraj, Theodore Arapoglou, Corbin Quick, Yaowu Liu, Han Chen, Ryan Sun, Rounak Dey, Donna K. Arnett, Lawrence F. Bielak, Joshua C. Bis, Thomas W. Blackwell, John Blangero, Eric Boerwinkle, Donald W. Bowden, Jennifer A. Brody, Brian E. Cade, Matthew P. Conomos, Adolfo Correa, L. Adrienne Cupples, Joanne E. Curran, Paul S. de Vries, Ravindranath Duggirala, Barry I. Freedman, Harald H. H. Göring, Xiuqing Guo, Rita R. Kalyani, Charles Kooperberg, Brian G. Kral, Leslie A. Lange, Ani Manichaikul, Lisa W. Martin, Braxton D. Mitchell, May E. Montasser, Alanna C. Morrison, Take Naseri, Jeffrey R. O’Connell, Nicholette D. Palmer, Patricia A. Peyser, Bruce M. Psaty, Laura M. Raffield, Susan Redline, Alexander P. Reiner, Muagututi‘a Sefuiva Reupena, Kenneth M. Rice, Stephen S. Rich, Jennifer A. Smith, Kent D.
Taylor, Ramachandran S. Vasan, Daniel E. Weeks, James G. Wilson, Lisa R. Yanek, Wei Zhao, NHLBI Trans-Omics for Precision Medicine (TOPMed) Consortium, TOPMed Lipids Working Group, Jerome I. Rotter, Christen J. Willer, Pradeep Natarajan, Gina M. Peloso and Xihong Lin. (2021). **A framework for detecting noncoding rare variant associations of large-scale whole-genome sequencing studies**. _bioRxiv_. DOI: <a href="https://doi.org/10.1101/2021.11.05.467531">10.1101/2021.11.05.467531</a>.

Xihao Li*, Zilin Li*, Hufeng Zhou, Sheila M. Gaynor, Yaowu Liu, Han Chen, Ryan Sun, Rounak Dey, Donna K. Arnett, Stella Aslibekyan, Christie M. Ballantyne, Lawrence F. Bielak, John Blangero, Eric Boerwinkle, Donald W. Bowden, Jai G. Broome, Matthew P. Conomos, Adolfo Correa, L. Adrienne Cupples, Joanne E. Curran, Barry I. Freedman, Xiuqing Guo, George Hindy, Marguerite R. Irvin, Sharon L. R. Kardia, Sekar Kathiresan, Alyna T. Khan, Charles L. Kooperberg, Cathy C. Laurie, X. Shirley Liu, Michael C. Mahaney, Ani W. Manichaikul, Lisa W. Martin, Rasika A. Mathias, Stephen T. McGarvey, Braxton D. Mitchell, May E. Montasser, Jill E. Moore, Alanna C. Morrison, Jeffrey R. O'Connell, Nicholette D. Palmer, Akhil Pampana, Juan M. Peralta, Patricia A. Peyser, Bruce M. Psaty, Susan Redline, Kenneth M. Rice, Stephen S. Rich, Jennifer A. Smith, Hemant K. Tiwari, Michael Y. Tsai, Ramachandran S. Vasan, Fei Fei Wang, Daniel E. Weeks, Zhiping Weng, James G. Wilson, Lisa R. Yanek, NHLBI Trans-Omics for Precision Medicine (TOPMed) Consortium, TOPMed Lipids Working Group, Benjamin M. Neale, Shamil R. Sunyaev, Gonçalo R. Abecasis, Jerome I. Rotter, Cristen J. Willer, Gina M. Peloso, Pradeep Natarajan, & Xihong Lin. (2020). **Dynamic incorporation of multiple in silico functional annotations empowers rare variant association analysis of large whole-genome sequencing studies at scale**. _Nature Genetics_, _52_(9), 969-983. PMID: <a href="https://www.ncbi.nlm.nih.gov/pubmed/32839606">32839606</a>. PMCID: <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7483769/">PMC7483769</a>. DOI: <a href="https://doi.org/10.1038/s41588-020-0676-4">10.1038/s41588-020-0676-4</a>.
## License
This software is licensed under GPLv3.

![GPLv3](http://www.gnu.org/graphics/gplv3-127x51.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)
