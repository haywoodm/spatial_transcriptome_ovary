# Spatially resolved transcriptomic profiling of ovarian aging in mice

## Overview
Ovarian aging precedes that of any other mammalian organ and is the primary instigator of female age-related infertility. The biological mechanisms responsible for ovarian aging and subsequent reproductive decline remain unclear. Previous studies have been limited by their use of traditional, bulk RNA-sequencing approaches, which masks the dynamic and heterogeneous nature of the ovary. In this study, we spatially resolve the transcriptomic landscape of ovaries from young and aged outbred mice using molecular barcoding technology to localize sequencing data within tissue sections. In total, we defined eight main ovarian cell populations, which were then further broken down into unique sub-clusters by gene expression. All clusters were characterized by significant changes to their transcriptomic landscape between young and aged samples, indicating that early onset aging affects all ovarian compartments. Differential analysis revealed that many pathways implicated in other age-related diseases, including Sirtuins Signaling and mTOR Signaling, are similarly dysregulated during ovarian aging. Further analysis of sub-cluster populations elucidated, for the first time, separate transcriptomes for distinctive granulosa cell sub-populations found in young and aged mice. Cell-type specific changes in aged granulosa cells were associated with the disruption of inhibins, activins, and gap junction activity. Additionally, four oocyte sub-cluster populations were uncovered, but only three were represented by aged ovaries. Overall, this study provides analysis of mammalian ovarian aging using spatial transcriptomics to achieve deeper understanding of the localization and cell-population-specific mechanisms underlying age-related fertility decline.

## Data files
Raw sequencing data, images, and processed data have been deposited in the GEO repository under accession GSE188257.

## R Scripts
Code to reproduce the analyses and figures from the manuscript are presented here.
- 1_preprocessing
   - partition slides by sample
   - add metadata
   - process feature matrix in Seurat
   - remove non-tissue data points
   - integrate slides
- 2_global_analysis
   - Figure 1b, 1c, 1d, 1e
   - Figure 2a
- 3_new_cluster_ids
   - Figure 3
   - Figure 2b
   - Figure 4
   - Table 1
   - Supplementary Table 1
- 4_subclusters
   - Figure 7 
   - Table 3
   - Supplementary Figure 3
   - Supplementary Figure 4
   - Supplementary Figure 5
- 5_oocytes_subcluster
   - Figure 8
