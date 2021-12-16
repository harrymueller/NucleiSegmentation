### Adapted BGI script

options(future.globals.maxSize = 100000000 * 1024^2,scipen = 20000000)
set.seed(6)

source("rscripts/functions/gem2rds.R")

library(argparser)
args <- arg_parser("Loading BGI spatial data(.gem)")
args <- add_argument(args, "--infile", help = "input .gem file")
args <- add_argument(args, "--outdir", help = "output directory")
args <- add_argument(args, "--binsize", help = "", default = 50)
args <- add_argument(args, "--proname", help = "", default = "None")
args <- add_argument(args, "--image", help = "Whether to generate the 'fake' image", default = NULL)
argv <- parse_args(args)

#if (!dir.exists(outdir)) dir.create(outdir)
#  outdir <- file.path(argv$outdir, argv$project_name)

if (argv$proname == "None") {
    argv$proname <- stringr::str_replace(tail(unlist(strsplit(argv$infile,"/")),1), ".gem", "")
    #project_name <- gsub(".bin.*","",project_name)
}

project_name <- paste0(argv$proname, "_bin", argv$binsize)
image = argv$image

seurat_spatialObj <- LoadBGI_Spatial(argv$infile, outdir = argv$outdir, bin_size = argv$binsize,  pro_name = project_name, UMI_GreyScale_Image = image)
