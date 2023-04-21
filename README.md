# Historical DNA reveals climate adaptation in an endangered songbird

This repository contains code to perform analyses described in the following paper:

Turbek SP, Bossu C, Rayne C, Gruppi C, Kus BE, Whitfield M, Smith TB, Paxton EH, Bay RA, and Ruegg KC (2023). Historical DNA reveals climate adaptation in an endangered songbird. Nature Climate Change.

## Abstract

To cope with climate change, species may shift their distributions or adapt in situ to changing environmental conditions. However, clear examples of genetic changes via adaptation are limited. We explore evolutionary responses to climate change in the endangered southwestern willow flycatcher (*Empidonax trailli extimus*) through whole-genome comparisons between historical specimens, collected from 1888-1909 near San Diego, CA, and contemporary individuals from across the breeding range. Genomic analyses revealed that introgression into San Diego increased adaptive potential over time and shifted genome-wide population structure towards that of neighboring populations. In contrast, loci linked to climate (dew point temperature and precipitation) shifted away from neighboring populations and in a direction consistent with adaptation to climate change in southern CA. This research highlights the role of admixture in facilitating adaptive shifts through its impact on genome-wide genetic variation and represents one of the few studies to document climate adaptation in a wild population.

[Turbek_Figure1.pdf](https://github.com/sturbek/WIFL-Climate-Adaptation/files/11298117/Turbek_Figure1.pdf)

Figure 1. Sampling Design and Patterns of Genomic Differentiation from Whole-Genome Sequencing Data (A) Sampling locations and subspecies boundaries across the willow flycatcher breeding range. Contemporary samples of the four subspecies are shown as squares, circles, diamonds, and triangles, while upside-down triangles indicate historical samples. (B) Genome-wide principal component analysis (PCA) (n = 238 individuals genotyped at 128,775 loci). Points are colored by sampling location and correspond to the map in (A). (C) Individual admixture proportions of contemporary and historical samples for K=4 (n = 238 individuals genotyped at 128,775 loci). Colored bars along the x axis indicate sampling locations and correspond to the map in (A).

## Code

Gradient_Forest_Pop_AllFreq_Extimus.R – Code to replicate the gradient forest analysis.
LFMM_Analysis_Indiv_Genotypes.R – Code to identify candidate climate-linked loci with LFMM.
RDA_Analysis_Indiv_Genotypes.R – Code to identify candidate climate-linked loci with RDA.

## Datasets

The data generated for this study are available in the Dryad Digital Repository (http://datadryad.org) and NCBI Sequence Read Archive (BioProject-PRJNA957938).
