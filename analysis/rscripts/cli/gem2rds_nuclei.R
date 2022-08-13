###
 # Converts segmented data into a Visium object
 # Adapted from a BGI script
###
library(argparser)
library(Matrix)

args <- arg_parser("Loading BGI spatial data(.gem)")
args <- add_argument(args, "--indir", help = "input .gem file") # assumes it has required files
args <- add_argument(args, "--outdir", help = "output directory")
argv <- parse_args(args)

assay = "Spatial"
slice = "stomics"

# parameters
dir = argv$indir # "/mnt/data/R_analysis/nuclei_seg"
filename = paste0(dir, "/nuclei.gem")
coords = paste0(dir, "/nuclei_coords.tsv")
spots = paste0(dir, "/nuclei_spots.tsv")
outdir = argv$outdir # "/mnt/data/R_analysis/nuclei_seg"

# read data
data <- data.table::fread(file = filename)
colnames(data) <- c("geneID", "nuclei_id", "UMICount")

coords = data.table::fread(file = coords)
spots = data.table::fread(file = spots)
  
# add mean spot coords to dat
data = merge(data, coords, by = "nuclei_id")

# converting to sparse matrix
gen_indices = match(data$geneID, unique(data$geneID))
nc_id_indices = match(data$nuclei_id, unique(data$nuclei_id))
mat <- Matrix::sparseMatrix(i = gen_indices, j = nc_id_indices, x = data$UMICount)

rownames(mat) = unique(data$geneID)
colnames(mat) = unique(data$nuclei_id)

# create seurat obj using matrix
object <- Seurat::CreateSeuratObject(mat, project = assay, assay = assay)

generate_image <- function(mat) {
  umiCountsPerNuclei = colSums(mat)
  scale_grey <- umiCountsPerNuclei / quantile(umiCountsPerNuclei, probs = 0.95)
  scale_grey[scale_grey > 1] <- 1
  
  tissue_df = data.frame(val = scale_grey, nuclei_id = colnames(mat))
  tissue_df = merge(tissue_df, spots, by="nuclei_id")
  tissue_lowres_image <- as.matrix(Matrix::sparseMatrix(
    i = tissue_df$y,
    j = tissue_df$x,
    x = tissue_df$val
  ))
  
  tissue_lowres_image_r <- raster::raster(tissue_lowres_image)
  tissue_lowres_image_r <- raster::writeRaster(tissue_lowres_image_r, file.path(outdir, "UMIGreyScaleFakeIMG.tiff"), overwrite=T)
  return(tissue_lowres_image)
}
im = generate_image(mat)

# define variables required for visium obj
tissue.positions <- data.frame(
  row.names = coords$nuclei_id,
  tissue = 1,
  row = coords$y, col = coords$x,
  imagerow = coords$y, imagecol = coords$x
)
scale.factors = list(
  fiducial_diameter_fullres = 1,
  tissue_hires_scalef = 1,
  tissue_lowres_scalef = 1
)
unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
spot.radius <- unnormalized.radius / max(dim(im)) # TODO possibly increases radius

# create visium object
spatialObj = new(Class = "VisiumV1",
          image = im,
          scale.factors = Seurat::scalefactors(
            spot = scale.factors$tissue_hires_scalef,
            fiducial = scale.factors$fiducial_diameter_fullres,
            hires = scale.factors$tissue_hires_scalef,
            lowres = scale.factors$tissue_lowres_scalef
          ),
          coordinates = tissue.positions,
          spot.radius = spot.radius
)

# add visium obj to seurat obj
spatialObj <- spatialObj[Seurat::Cells(object)]
Seurat::DefaultAssay(spatialObj) <- assay
object[[slice]] <- spatialObj

object <- subset(object, subset = nCount_Spatial > 0)

# save spot maps with object - used for plotting
colnames(coords) = c("nuclei_id", "cx", "cy")
data = list("seurat" = object, "spot_mappings" = merge(coords, spots, by = "nuclei_id"))

# save file
saveRDS(data, file = file.path(outdir, "spatialObj.rds"))

print("Completed successfully")