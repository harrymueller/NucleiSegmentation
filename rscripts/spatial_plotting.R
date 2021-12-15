# LIBRARIES
library(Seurat)
library(ggplot2)
library(patchwork)
library(raster)
library(RColorBrewer)

# params
library(argparser)
args <- arg_parser("Plotting nCount, nFeature, Malat1 and Neat1 violin and spatial plots.")
args <- add_argument(args, "--binsize", help = "Bin Size")
args <- add_argument(args, "--id", help = "Tongue ID")
argv <- parse_args(args)

BIN_SIZE = argv$binsize
TONGUE_ID = argv$id

# other consts
INPUT_DIR = "/mnt/data/tongue/discovery/gemRDS"
OUTPUT_DIR = sprintf("/mnt/data/tongue/discovery/count_feature_plots/%s_bin%s", TONGUE_ID, BIN_SIZE)
source("/mnt/data/tongue/scripts/rscripts/accurate_plot.R")

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

# takes SpatialPlot data from seurat and saves it to a file


# read in RDS
INPUT = sprintf("%s/%s_bin%s_spatialObj.rds", INPUT_DIR, TONGUE_ID, BIN_SIZE)
obj = readRDS(INPUT)
gc()

# plots
colours = c("white",rev(brewer.pal(n = 11, name = "Spectral")))

for (t in c("nCount_Spatial", "nFeature_Spatial", "Malat1", "Neat1")) {
  p = VlnPlot(obj, features = t, pt.size = 0.05) + NoLegend()
  isGene = t %in% c("Malat1", "Neat1")
  accurate_plot(SpatialPlot(obj, features = "nCount_Spatial")$data, 
                filename = paste0(OUTPUT_DIR, "/nCount_Spatial.png"), 
                legend_name = paste0("nCounts", isGene ? "" : " (q0.99999)",
                adjust = isGene ? 1 : 0.9999,
                custom_colours = colours,
                left_plot = p)
  gc()
} 


