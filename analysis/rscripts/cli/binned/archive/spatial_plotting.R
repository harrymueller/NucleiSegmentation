# LIBRARIES
library(Seurat)
library(ggplot2)
library(patchwork)
library(raster)
library(RColorBrewer)
library(viridis)

# params
library(argparser)
args <- arg_parser("Plotting nCount, nFeature, Malat1 and Neat1 violin and spatial plots.")
args <- add_argument(args, "--input", help = "Input")
args <- add_argument(args, "--output", help = "Output dir")
argv <- parse_args(args)

INPUT = argv$input
OUTPUT_DIR = argv$output

# INPUT OUTPUT DIR
GENES = c("Krt36", "Krt84", "Krt35", "Krt6b", "Krt6a", "Dmd", "Neb", "Ttn", "Rbfox1", "C1qtnf7", "Dclk1", "Mmp16", "Mecom", "Plcb1", "Flt1", "Fetub", "Ecm1", "Fmo2", "Pappa", "Hs3st5", "St3gal4", "Mki67", "Cenpf", "Cenpe", "Frem2", "Sema3c", "Pcdh7", "Unc5c", "Cpm", "Hhip", "Hephl1") 

source("/mnt/data/scripts/analysis/rscripts/functions/accurate_plot.R")
  
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

# read in RDS
obj = readRDS(INPUT)

# plots
#colours = c("white",rev(brewer.pal(n = 11, name = "Spectral")))
colours = rev(brewer.pal(n = 11, name = "Spectral"))
i = 0
#for (t in c("nCount_Spatial", "nFeature_Spatial", "Malat1", "Neat1")) {
for (t in GENES) {
  if (BIN_SIZE > 11) {
    #left_plot = VlnPlot(obj, features = t, pt.size = 0.05) + NoLegend()
    isGene = t %in% GENES
    accurate_plot(SpatialPlot(obj, features = t)$data, 
                  filename = paste0(OUTPUT_DIR, "/", t, ".png"), 
                  legend_name = paste0(t, ifelse(!isGene && BIN_SIZE <= 10, "\n(q0.99999)", "")),
                  adjust = ifelse(!isGene && BIN_SIZE <= 10, 0.9999, 1),
                  custom_colours = colours, #left_plot = left_plot,
                  dpi = 750, # increasing DPI increases size of scales etc.
                  minres = 1500)
  } else {
    # plot vln and spatial plot separately for bin <= 2
    #left_plot = VlnPlot(obj, features = t, pt.size = 0.05) + NoLegend()
    #ggsave(paste0(OUTPUT_DIR, "/", t, "_vln.png"), 
    #       left_plot,
    #       width = 4)
    isGene = t %in% GENES
    accurate_plot(SpatialPlot(obj, features = t)$data, 
                  filename = paste0(OUTPUT_DIR, "/", t, "_spatial.png"), 
                  legend_name = paste0(t, ifelse(!isGene && BIN_SIZE <= 10, "\n(q0.99999)", "")),
                  adjust = ifelse(!isGene && BIN_SIZE <= 10, 0.9999, 1),
                  custom_colours = colours,
                  dpi = 750, #2000,
                  minres = 1500, #1,
                  crop = TRUE,
                  black_background = TRUE)
  }
  gc()
} 


