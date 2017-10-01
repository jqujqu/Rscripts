
library("optparse")

opt.list = list(
    make_option(c("-s", "--srainfo"),
                type    = "character",
                default = "",
                help    = "sra run information table (output csv from sraInfo.R)"),
    make_option(c("-f", "--field"),
                type    = "character",
                default = "sample",
                help    = paste("column name of the field in the SRAINFO table to rename,",
                                "Default is 'sample'.")),
    make_option(c("-o", "--outdir"),
                type    = "character",
                default = getwd(),
                help    = paste("Path to create subdir structure.",
                                "Default is CWD."))
)

# Read options from command line
parser <- OptionParser(usage = "%prog [options] field_rename_table \n#(rename convention: <Species>_<CellType>_R<rep>_L<rep>)",
                       option_list = opt.list)
arguments  <- parse_args(parser, positional_arguments = c(1, Inf))
opt        <- arguments$options
renamefile <- arguments$args[1]

x <- read.csv(opt$srainfo, as.is=T)
y <- read.table(renamefile, as.is=T)

x$rename <- ""
for(i in 1:nrow(y)) {
 x$rename[which(x[, opt$field]==y$V1[i])] <- y$V2[i]
}

for (i in unique(x$rename)){
  parts <- unlist(strsplit(i, "_"))
  sample <- paste(parts[1],parts[2], sep="_")
  biorep <- paste(parts[1], parts[2], parts[3], sep="_")
  techrep <- i
  path <- paste(opt$outdir, sample, biorep, techrep, sep="/")
  if (!dir.exists(path)) {
    message(paste("Creating path ", path))
    dir.create(path, recursive=T)
  }
  message(paste("Writing related run information to ", paste(path, "runinfo.csv", sep="/")))
  write.csv(x[x$rename == i, ], file=paste(path, "runinfo.csv", sep="/"))
  write.table(x$ftp[x$rename == i], quote=F, row.names=F, col.names=F,
              file=paste(path, "sralinks.txt", sep="/"))
}
