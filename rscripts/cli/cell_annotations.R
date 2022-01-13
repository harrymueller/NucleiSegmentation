##############################
# LIBS + PARAMS
##############################
library(Seurat)
library(SingleR)
library(ggplot2)

RefData = celldex::MouseRNAseqData()

DIR = "/mnt/data"

# to do multiple runs in 1 go
DATA = data.frame(
  "id"          = c("tongue-5"),
  "binsize"     = c(50),
  "res"         = c(0.4),
  "norm_method" = c("SCTransform"),
  "diameter"    = c(0)
)

##############################
# FUNCTIONS
##############################
singleR_annotations <- function (bUseMainLabels) {
    
}

##############################
# START
##############################
for (it in seq(1, length(rownames(DATA)))) {
    current = DATA[it,]
    print("##############################")
    print(paste0(
        "Starting ", current$id, "_bin", current$binsize, 
        ifelse(current$diameter == 0, "", paste0("_subset",current$diameter)), " ..."
    ))
    print("##############################")

    # create filenames
    if (current$diameter == 0) { # not a subset
        filename = sprintf("%s_bin%s_red.Rds", current$id, current$binsize)
        folder_name = sprintf("%s_bin%s", current$id, current$binsize)
    } else {
        filename = sprintf("%s_bin%s_subset%s_red.Rds", current$id, current$binsize, current$diameter)
        folder_name = sprintf("%s_bin%s_subset%s", current$id, current$binsize, current$diameter)
    }

    # output dir
    output_dir = sprintf("%s/cell_annotations/%s", DIR, folder_name)
    if (!dir.exists(output_dir)) dir.create(output_dir)

    # read data
    input = sprintf(
        "%s/%s/%s", DIR, ifelse(current$norm_method == "SCTransform", "scDimReducedRDS", "dimReducedRDS"), filename
    )
    obj = readRDS(input)
    
    assay_to_use = ifelse(current$norm_method == "SCTransform", "SCT", "Spatial")

    



    # done
    print("##############################")
    print(paste0(
        "Finished ", current$id, "_bin", current$binsize, 
        ifelse(current$diameter == 0, "", paste0("_subset",current$diameter))
    ))
    print("##############################")

    gc()
}