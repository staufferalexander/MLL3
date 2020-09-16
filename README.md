# MLL3
This repository hosts the code used in my research on how the loss of MLL3 affects ER-alpha regulation in ER+ breast cancer. Research performed at Vanderbilt University, under Thomas Stricker, MD PhD.

1. GERP analysis

1. DEMETER dependence tool - https://depmap.org/rnai/genedeps?gene=KMT2C

1. RNA-seq
   1. Tophat
      * Version 2.0.13
   1. Cufflinks
      * Version 2.2.1
   1. Cuffnorm  
      * Version 2.2.1
   1. Differential Gene Expression Scripts
      1. TCGA
      1. ZR751
   1. WebGestalt
      * WebGestalt 2019 Version
      * Parameters
         * Minimum Number of Genes for a Category = 3
         * Maximum Number of Genes for a Category = 2000
         * Multiple Test Adjustment = BH
         * Functional Database = others
            * Uploaded MSigDB Release 7.0
   1. iRegulon
      * Cytoscape Version 3.7.1 
1. ChIP-seq
   1. BWA Alignment
      * Version 0.7.5a-r405
   1. Peak Calling via SPP according to ENCODE best practices
      * SPP version 1.15.5
      * ENCODE 3 Pipeline v1 https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#
      * IDR version 2.0.3
   1. Enrichment Analysis
      1. GREAT
         * Version 3.0.0
      1. WebGestalt
         * See Info Above
   1. Motif Analysis
      1. MEME-Suite
         * Version 4.11.2
   1. Peak Assignment
      1. Bedtools
         * Version 2.26.0
      1. Assessment of permutation-based assignment
1. Proliferation Assays
   1. GRMetrics
      * Version 1.10.0
      * R Version 3.6.1
      * Bioconducter Version 3.9
1. Accession Numbers
    1. RNA-seq data
    1. ChIP-seq data
