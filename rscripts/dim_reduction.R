options(future.globals.maxSize = 100000000 * 1024^2,scipen = 20000000)
set.seed(6)

library(argparser)
args <- arg_parser("Loading BGI spatial data(.gem)")
args <- add_argument(args, "--binsize", help = "binsize")
args <- add_argument(args, "--id", help="id")
argv <- parse_args(args)


OUTPUT = sprintf("/data/tongue/dimReducedRDS/%s_bin%s_red.Rds", TONGUE_ID, BIN_SIZE)

DIR = "/mnt/data/tongue_STOmics/discovery/gemRDS"
INPUT = sprintf("%s/%s_bin%s_spatialObj.rds", DIR, TONGUE_ID, BIN_SIZE)
obj = readRDS(INPUT)

# method of normalisation
obj <- SCTransform(obj, assay = "Spatial", verbose = FALSE)

# dimension reduction
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)

# save new obj
saveRDS(obj, OUTPUT)