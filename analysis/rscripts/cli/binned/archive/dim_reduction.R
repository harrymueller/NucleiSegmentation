options(future.globals.maxSize = 100000000 * 1024^2,scipen = 20000000)
set.seed(6)

library(argparser)
library(Seurat)

args <- arg_parser("Loading BGI spatial data(.gem)")
args <- add_argument(args, "--binsize", help = "binsize")
args <- add_argument(args, "--id", help="id")
args <- add_argument(args, "--method", help="method of normalisation, SCT || LN")
args <- add_argument(args, "--diameter", help="if subsetting, supply the diameter", default=NULL)

argv <- parse_args(args)

TONGUE_ID = argv$id
BIN_SIZE = argv$binsize
METHOD = argv$method
DIAMETER = argv$diameter

METHOD_NAME = ifelse(METHOD == "SCT", "SCTransform", "NormalizeData")

if (DIAMETER == 0) {
  INPUT_DIR = "/mnt/data/R_analysis/gemRDS"
  INPUT = sprintf("%s/%s_bin%s_spatialObj.rds", INPUT_DIR, TONGUE_ID, BIN_SIZE)
} else {
  INPUT_DIR = "/mnt/data/R_analysis/subsets"
  INPUT = sprintf("%s/%s_bin%s_subset%s.rds", INPUT_DIR, TONGUE_ID, BIN_SIZE, DIAMETER)
}

obj = readRDS(INPUT)    

if (METHOD == "SCT") {
    if (DIAMETER == 0)
        OUTPUT = sprintf("/mnt/data/R_analysis/dimReduction/%s/%s_bin%s_red.Rds", METHOD_NAME, TONGUE_ID, BIN_SIZE)
    else
        OUTPUT = sprintf("/mnt/data/R_analysis/dimReduction/%s/%s_bin%s_subset%s_red.Rds", METHOD_NAME, TONGUE_ID, BIN_SIZE, DIAMETER)
    
    obj <- SCTransform(obj, assay = "Spatial", verbose = FALSE)
    obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
} else if (METHOD == "LN") {
    if (DIAMETER == 0) {
        OUTPUT = sprintf("/mnt/data/R_analysis/dimReducedRDS/%s_bin%s_red.Rds", TONGUE_ID, BIN_SIZE)
    } else {
        OUTPUT = sprintf("/mnt/data/R_analysis/dimReducedRDS/%s_bin%s_subset%s_red.Rds", TONGUE_ID, BIN_SIZE, DIAMETER)
    }
    
    obj <- NormalizeData(obj, assay = "Spatial", verbose = FALSE)

    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, assay = "Spatial", verbose = FALSE)
} else {
    sprintf("Unrecognised parameter for --method: %s. Use SCT or LN", METHOD)
    exit()
}

# dimension reduction
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)

# save new obj
saveRDS(obj, OUTPUT)