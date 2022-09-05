library(argparser)
args <- arg_parser("Loading BGI spatial data(.gem)")
args <- add_argument(args, "--infile", help = "input .gem file")
args <- add_argument(args, "--outdir", help = "output directory")
args <- add_argument(args, "--binsize", help = "", default = 50)
args <- add_argument(args, "--proname", help = "", default = "None")
args <- add_argument(args, "--image", help = "Whether to generate the 'fake' image", default = NULL)
argv <- parse_args(args)

filename = argv$infile
project_name <- paste0(argv$proname, "_bin", argv$binsize)
outdir = argv$outdir
bin_size = argv$binsize
pro_name = project_name


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

print("binned")
dat$bin_ID <- max(dat$x) * (dat$y - 1) + dat$x
bin.coor <- dat[, sum(V1), by = .(x, y)]

scale_grey <- bin.coor$V1 / quantile(bin.coor$V1, probs = 0.95)
scale_grey[scale_grey > 1] <- 1
tissue_lowres_image <- as.matrix(Matrix::sparseMatrix(
    i = bin.coor$y,
    j = bin.coor$x,
    x = scale_grey
))
tissue_lowres_image_r <- raster::raster(tissue_lowres_image)
tissue_lowres_image_r <- raster::writeRaster(tissue_lowres_image_r, file.path(outdir, paste0(pro_name, "_UMIGreyScaleFakeIMG.tiff")), overwrite=T)
