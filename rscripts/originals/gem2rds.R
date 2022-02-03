options(future.globals.maxSize = 100000000 * 1024^2,scipen = 20000000)
set.seed(6)

library(argparser)
args <- arg_parser("Loading BGI spatial data(.gem)")
args <- add_argument(args, "--infile", help = "input .gem file")
args <- add_argument(args, "--outdir", help = "output directory")
args <- add_argument(args, "--binsize", help = "", default = 50)
args <- add_argument(args, "--proname", help = "", default = "None")
argv <- parse_args(args)

infile <- argv$infile
outdir <- argv$outdir
binsize <- argv$binsize
project_name <- argv$proname



if (!dir.exists(outdir)) dir.create(outdir)
outdir <- file.path(outdir, project_name)
if (!dir.exists(outdir)) dir.create(outdir)
saveRDS(argv, file.path(outdir, "arguments.rds"))

if (project_name == "None") {
    project_name <- tail(unlist(strsplit(infile,"/")),1)
    project_name <- gsub(".bin.*","",project_name)
}

project_name <- paste0(project_name, "_bin", binsize)

LoadBGI_Spatial <- function(filename,
                            outdir = getwd(),
                            bin_size = 50,
                            save_rds = T,
                            pro_name = "Spatial",
                            UMI_GreyScale_Image = T,
                            assay = "Spatial",
                            slice = "slice1") {
    
    dat <- data.table::fread(file = filename)
    colnames(dat)[1:4] <- c("geneID", "x", "y", "UMICount")
    rowname_style <- "BIN"

    dat$x <- trunc((dat$x - min(dat$x)) / bin_size + 1)
    dat$y <- trunc((dat$y - min(dat$y)) / bin_size + 1)


    if ("MIDCounts" %in% colnames(dat)) {
        dat <- dat[, sum(MIDCounts), by = .(geneID, x, y)]
    } else if ("MIDCount" %in% colnames(dat)) {
	    dat <- dat[, sum(MIDCount), by = .(geneID, x, y)]	
    } else {
        dat <- dat[, sum(UMICount), by = .(geneID, x, y)]
    }
    dat$bin_ID <- max(dat$x) * (dat$y - 1) + dat$x
    bin.coor <- dat[, sum(V1), by = .(x, y)]
    if (UMI_GreyScale_Image) {
        scale_grey <- bin.coor$V1 / quantile(bin.coor$V1, probs = 0.95)
        scale_grey[scale_grey > 1] <- 1
        tissue_lowres_image <- as.matrix(Matrix::sparseMatrix(
        i = bin.coor$y,
        j = bin.coor$x,
        x = scale_grey
        ))
        tissue_lowres_image_r <- raster::raster(tissue_lowres_image)
        tissue_lowres_image_r <- raster::writeRaster(tissue_lowres_image_r, file.path(outdir, paste0(pro_name, "_UMIGreyScaleFakeIMG.tiff")), overwrite=T)
    } else tissue_lowres_image <- matrix(0, max(bin.coor$y), max(bin.coor$x))
    
    geneID <- seq(length(unique(dat$geneID)))
    hash.G <- data.frame(row.names = unique(dat$geneID), values = geneID)
    gen <- hash.G[dat$geneID, "values"]
    bin_ID <- unique(dat$bin_ID)
    hash.B <- data.frame(row.names = sprintf("%d", bin_ID), values = bin_ID)
    bin <- hash.B[sprintf("%d", dat$bin_ID), "values"]
    cnt <- dat$V1
    mat <- Matrix::sparseMatrix(i = gen, j = bin, x = cnt)
    rownames(mat) <- rownames(hash.G)
    colnames(mat) <- paste(rowname_style, sprintf("%d", seq(max(hash.B[, "values"]))), sep = ".")
    object <- Seurat::CreateSeuratObject(mat, project = assay, assay = assay)
    rm(dat)
    gc()

    generate_spatialObj <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE) {
        if (filter.matrix) {
            tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
        }
        unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
        spot.radius <- unnormalized.radius / max(dim(image))
        return(new(Class = "VisiumV1",
                    image = image,
                    scale.factors = Seurat::scalefactors(
                        spot = scale.factors$tissue_hires_scalef,
                        fiducial = scale.factors$fiducial_diameter_fullres,
                        hires = scale.factors$tissue_hires_scalef,
                        lowres = scale.factors$tissue_lowres_scalef
                    ),
                    coordinates = tissue.positions,
                    spot.radius = spot.radius
                ))
    }

    tissue_positions_list <- data.frame(
        row.names = paste(rowname_style, rownames(hash.B), sep = "."),
        tissue = 1,
        row = bin.coor$y, col = bin.coor$x,
        imagerow = bin.coor$y, imagecol = bin.coor$x
    )
    scalefactors_json <- rjson::toJSON(list(
        fiducial_diameter_fullres = 1,
        tissue_hires_scalef = 1,
        tissue_lowres_scalef = 1
    ))
    spatialObj <- generate_spatialObj(
        image = tissue_lowres_image,
        scale.factors = rjson::fromJSON(scalefactors_json),
        tissue.positions = tissue_positions_list
    )

    spatialObj <- spatialObj[Seurat::Cells(object)]
    Seurat::DefaultAssay(spatialObj) <- assay
    object[[slice]] <- spatialObj

    object <- subset(object, subset = nCount_Spatial > 0)

    if (save_rds) {
        saveRDS(object, file = file.path(outdir, paste0(pro_name, "_spatialObj.rds")))
    }
    return(object)
}
seurat_spatialObj <- LoadBGI_Spatial(infile, outdir = outdir, bin_size = binsize,  pro_name = project_name)
