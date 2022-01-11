library(Seurat)

# params
library(argparser)
args <- arg_parser("Plotting nCount, nFeature, Malat1 and Neat1 violin and spatial plots.")
args <- add_argument(args, "--binsize", help = "Bin Size")
args <- add_argument(args, "--id", help = "Tongue ID")
args <- add_argument(args, "--diameter", help = "Diameter - num bins from outer contour")
argv <- parse_args(args)

BIN_SIZE = argv$binsize
TONGUE_ID = argv$id
DIAMETER = argv$diameter

INPUT_DIR = "/mnt/data/gemRDS"
SUBSET_FILE = sprintf("/mnt/data/subsets/%s_bin%s_subset%s.tsv", TONGUE_ID, BIN_SIZE, DIAMETER)
OUTPUT_FILE = sprintf("/mnt/data/subsets/%s_bin%s_subset%s.rds", TONGUE_ID, BIN_SIZE, DIAMETER)

obj = readRDS(sprintf("%s/%s_bin%s_spatialObj.rds", INPUT_DIR, TONGUE_ID, BIN_SIZE))

tab = read.csv(SUBSET_FILE, sep="\t")
# fixed from 0-indexed to 1-indexed
tab[[1]] = tab[[1]] + 1
tab[[2]] = tab[[2]] + 1 

# create df containing "y x" with rownames bin IDs 
coords = obj@images$slice1@coordinates
coords = data.frame("ids" = rownames(coords),
    "cat" = paste(coords$row, coords$col))

# combine y and x for new coords
tab = paste(tab$y_coord, tab$x_coord)

# create a list of ids to keep
subset_ids = coords$ids[coords$cat %in% tab]

# add bin_id and subset (T || F) to metadata
obj@meta.data[["bin_id"]] = names(Idents(obj))
obj@meta.data[["subset"]] = names(Idents(obj)) %in% subset_ids

# subset
print("Subsetting...")
Idents(obj) = "subset"
new_obj = subset(obj, idents = TRUE)
Idents(new_obj) = "bin_id"

new_obj$subset = NULL

# save
print("Saving...")
saveRDS(new_obj, file = OUTPUT_FILE)
