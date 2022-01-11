options(future.globals.maxSize = 100000000 * 1024^2,scipen = 20000000)
set.seed(6)

library(argparser)
library(Seurat)

args <- arg_parser("Loading BGI spatial data(.gem)")
args <- add_argument(args, "--binsize", help = "binsize")
args <- add_argument(args, "--id", help="id")
args <- add_argument(args, "--method", help="method of normalisation, SCT || LN")
argv <- parse_args(args)

TONGUE_ID = argv$id
BIN_SIZE = argv$binsize
METHOD = argv$method

DIR = "/mnt/data/gemRDS"
INPUT = sprintf("%s/%s_bin%s_spatialObj.rds", DIR, TONGUE_ID, BIN_SIZE)
obj = readRDS(INPUT)    

if (METHOD == "SCT") {
    OUTPUT = sprintf("/mnt/data/scDimReducedRDS/%s_bin%s_red.Rds", TONGUE_ID, BIN_SIZE)
    obj <- SCTransform(obj, assay = "Spatial", verbose = FALSE)

    obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
} else if (METHOD == "LN") {
    OUTPUT = sprintf("/mnt/data/dimReducedRDS/%s_bin%s_red.Rds", TONGUE_ID, BIN_SIZE)
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
obj <- FindClusters(obj, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)

# save new obj
saveRDS(obj, OUTPUT)