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
if (T) {
  INPUT_DIR = "/mnt/data/gemRDS"
  OUTPUT_DIR = sprintf("/mnt/data/count_feature_plots/%s_bin%s", TONGUE_ID, BIN_SIZE)
  source("/mnt/data/scripts/rscripts/accurate_plot.R")
} else {
  INPUT_DIR = "/mnt/data/tongue_STOmics/discovery/gemRDS"
  OUTPUT_DIR = sprintf("/mnt/data/tongue_STOmics/discovery/count_feature_plots/%s_bin%s", TONGUE_ID, BIN_SIZE)
  source("/mnt/local/scripts/rscripts/accurate_plot.R")
}
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

# read in RDS
INPUT = sprintf("%s/%s_bin%s_spatialObj.rds", INPUT_DIR, TONGUE_ID, BIN_SIZE)
obj = readRDS(INPUT)

# plots
#colours = c("white",rev(brewer.pal(n = 11, name = "Spectral")))
colours = rev(brewer.pal(n = 11, name = "Spectral"))

for (t in c("nCount_Spatial", "nFeature_Spatial", "Malat1", "Neat1")) {
#for (t in c("nCount_Spatial")) {
  left_plot = VlnPlot(obj, features = t, pt.size = 0.05) + NoLegend()
  isGene = t %in% c("Malat1", "Neat1")
  accurate_plot(SpatialPlot(obj, features = t)$data, 
                filename = paste0(OUTPUT_DIR, "/", t, ".png"), 
                legend_name = paste0(t, ifelse(!isGene && BIN_SIZE <= 10, "\n(q0.99999)", "")),
                adjust = ifelse(!isGene && BIN_SIZE <= 10, 0.9999, 1),
                custom_colours = colours,
                left_plot = left_plot,
                dpi = 2000,
                minres = 1)
  gc()
} 


