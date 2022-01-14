##############################
# LIBS + PARAMS
##############################
library(Seurat)
library(SingleR)
library(ggplot2)

RefData = celldex::MouseRNAseqData()

# params
DIR = "/mnt/data"

# arg parsing
library(argparser)
args <- arg_parser("Cluster the seurat obj then plot a umap and spatial image")
args <- add_argument(args, "--id", help="TONGUE_ID")
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

# other required variables
METHOD_NAME = ifelse(METHOD == "SCT", "SCTransform", "NormalizeData")
ASSAY_TO_USE = ifelse(METHOD == "SCT", "SCT", "Spatial")

FOLDER_NAME = sprintf("%s_bin%s_subset%s_res%s", TONGUE_ID, BINSIZE, DIAMETER, RESOLUTION)
FILENAME = paste0(FOLDER_NAME, ".rds")
INPUT = paste0(DIR, "/umap_clusters/", METHOD_NAME, "/RDS/", FILENAME)

# output dir
OUTPUT_DIR = sprintf("%s/markers/%s/%s", DIR, METHOD_NAME, FOLDER_NAME)
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

##############################
# FUNCTIONS
##############################
singleR_annotations <- function (bUseMainLabels) {
    
}

##############################
# START
##############################
# read file
obj = readRDS(INPUT)

