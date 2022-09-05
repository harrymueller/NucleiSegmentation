# LIBRARIES
library(Seurat)
library(ggplot2)
library(patchwork)
library(raster)
library(RColorBrewer)

# params
library(argparser)
args <- arg_parser("Cluster the seurat obj then plot a umap and spatial image")
args <- add_argument(args, "--infile", help = "Input File")
args <- add_argument(args, "--outdir", help = "Output Dir")
args <- add_argument(args, "--resolution", help="seurat FindCluster resolution", default=0.5)
argv <- parse_args(args)

OUTPUT_DIR = argv$outdir
INPUT = argv$infile
RESOLUTION = as.double(argv$resolution)

source("../../functions/accurate_plot.R")

# read in RDS
data = readRDS(INPUT)
obj = data$seurat

# find clusters and plot
obj <- FindClusters(obj, resolution = RESOLUTION)
p1 <- DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 1)
p2 <- SpatialDimPlot(obj)

NAME = paste0("umap ", RESOLUTION)

accurate_plot(
    p2$data,
    filename = sprintf("%s/%s.png", OUTPUT_DIR, NAME),
    legend_name = stringr::str_replace_all(NAME, "_", "\n"),
    left_plot = p1,
    dpi = 500,
    minres = 1000,
    spot_mappings = data$spot_mappings
)

saveRDS(obj, sprintf("%s/%s.rds", OUTPUT_DIR, NAME))
