options(future.globals.maxSize = 100000000 * 1024^2,scipen = 20000000)
set.seed(6)

library(argparser)
library(Seurat)

args <- arg_parser("Loading BGI spatial data(.gem)")
args <- add_argument(args, "--infile", help = "Input File")
args <- add_argument(args, "--outdir", help = "Output Dir")
args <- add_argument(args, "--method", help="method of normalisation, SCT || LN")
argv <- parse_args(args)

OUTPUT_DIR = argv$outdir
INPUT = argv$infile
METHOD = argv$method

OUTPUT = paste0(OUTPUT_DIR, "/", "dim_red_", METHOD, ".rds")
METHOD_NAME = ifelse(METHOD == "SCT", "SCTransform", "NormalizeData")

data = readRDS(INPUT)    
obj = data$seurat

if (METHOD == "SCT") {
    obj <- SCTransform(obj, assay = "Spatial", verbose = FALSE)
    obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
} else if (METHOD == "LN") {
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