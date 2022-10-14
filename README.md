# Benchmarking Nuclei Segmentation Algorithms using Real and Synthetic Data for Improving STOmics Analyses

## Harrison Mueller
 
## Supervisors
#### Harry Perkins: *Systems Biology and Genomics Laboratory*
 - Elena Denisenko
 - Alistair Forrest
 
#### UWA: *Computer Science Department*
 - Amitava Datta
 - Max Ward
 
The code is this repository was used as part of my honours project. 

## Abstract
Spatial transcriptomics is the latest technology for studying transcriptomics, and is becoming an increasingly important technique for cellular biology, neuroscience, developmental biology, and cancer research. Several platforms have been developed to measure RNA expression as well as the spatial origin of a given transcript, with one of the most recent being STOmics (Spatial TranscriptOmics), formally Stereo-seq.

STOmics has a subcellular resolution, meaning that each spot where transcripts are recorded is smaller than a cell or nuclei. However, this has three main drawbacks: 1) spots are not representative of cells or nuclei; 2) the per-spot data is sparse; and 3) even small biological samples produce large amounts of data, which require significant amounts of computational power to analyse.

These issues can be solved by aggregating spots. Classically, a grid is used to do this, where any spots that fall within a square are aggregated into a single “cell”. A better way of doing this would be to use computer vision algorithms to segment nuclei based on a nuclei stained image of the sequenced tissue section. These segmentation maps can then be applied to STOmics data to obtain per-nuclei counts.
This solves the three problems as: 1) the per-nuclei counts contain expression levels of specific nuclei, 2) the per-nuclei data is not sparse, and 3) the size of the data is heavily reduced, requiring significantly less computational power to  analyse.

Various segmentation algorithms have been developed for use on similar image types and are therefore likely to segment ssDNA stains effectively. Here, we compare and contrast three such algorithms: two pre-trained neural networks, Cellpose and DeepCell; and a conventional segmentation algorithm, the watershed algorithm.

To compare and contrast these tools, the segmentation of ssDNA-stained images of mouse tongue samples were compared qualitatively. For a quantitative comparison, synthetic ssDNA stain images were generated such that the ground-truths were known, following which, performance metrics were computed. Furthermore, the alignment of the nuclei stained image to the transcriptional data was investigated via the distribution of pre-mRNA, which is only present in nuclei. The data has also been explored after aggregating spots using the classical binning approach.



