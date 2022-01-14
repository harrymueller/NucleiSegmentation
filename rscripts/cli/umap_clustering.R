# LIBRARIES
library(Seurat)
library(ggplot2)
library(patchwork)
library(raster)
library(RColorBrewer)

# params
library(argparser)
args <- arg_parser("Cluster the seurat obj then plot a umap and spatial image")
args <- add_argument(args, "--id", help="id")
args <- add_argument(args, "--binsize", help = "binsize")
args <- add_argument(args, "--method", help="method of normalisation, SCT || LN")
args <- add_argument(args, "--diameter", help="supply the diameter, 0 if not subsetting")
args <- add_argument(args, "--resolution", help="seurat FindCluster resolution", default=0.5)
argv <- parse_args(args)

TONGUE_ID   = argv$id
METHOD      = argv$method
BINSIZE     = as.integer(argv$binsize)
DIAMETER    = as.integer(argv$diameter)
RESOLUTION  = as.double(argv$resolution)

METHOD_NAME = ifelse(METHOD == "SCT", "SCTransform", "NormalizeData")
METHOD_FOLDER = ifelse(METHOD == "SCT", "scDimReducedRDS", "dimReducedRDS")

# other consts
INPUT_DIR = "/mnt/data/dimReduction"
NAME = sprintf("%s_bin%s_subset%s_res%s", TONGUE_ID, BINSIZE, DIAMETER, round(RESOLUTION, 1))
PLOTS_DIR = sprintf("/mnt/data/umap_clusters/%s/plots", METHOD_NAME)
RDS_DIR = sprintf("/mnt/data/umap_clusters/%s/RDS", METHOD_NAME)

source("/mnt/data/scripts/rscripts/functions/accurate_plot.R")

if (!dir.exists(PLOTS_DIR)) dir.create(PLOTS_DIR)
if (!dir.exists(RDS_DIR)) dir.create(RDS_DIR)

# read in RDS
if (DIAMETER == 0) {
    INPUT = sprintf("%s/%s/%s_bin%s_red.Rds", INPUT_DIR, METHOD_NAME, TONGUE_ID, BINSIZE)
} else {
    INPUT = sprintf("%s/%s/%s_bin%s_subset%s_red.Rds", INPUT_DIR, METHOD_NAME, TONGUE_ID, BINSIZE, DIAMETER)
}
obj = readRDS(INPUT)

# find clusters and plot
obj <- FindClusters(obj, resolution = RESOLUTION)
p1 <- DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 1)
p2 <- SpatialDimPlot(obj)

accurate_plot(
    p2$data,
    filename = sprintf("%s/%s.png", PLOTS_DIR, NAME),
    legend_name = stringr::str_replace_all(NAME, "_", "\n"),
    left_plot = p1,
    dpi = 500,
    minres = 1000
)

saveRDS(obj, sprintf("%s/%s.rds", RDS_DIR, NAME))
