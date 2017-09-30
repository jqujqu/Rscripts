library("optparse")

opt.list = list(
    make_option(c("-d", "--database"),
                type    = "character",
                default = "./SRAmetadb.sqlite",
                help    = paste("Path to SRA database file, if does not exist",
                                "database will be download to this location.",
                                "Default is './SRAmetadb.sqlite'.")),
    make_option(c("-m", "--geometadb"),
                type    = "character",
                default = "./GEOmetadb.sqlite",
                help    = paste("Path to GEOmetadb database file, if does not exist",
                                "database will be download to this location.",
                                "Default is './GEOmetadb.sqlite'.")),
    make_option(c("-o", "--outdir"),
                type    = "character",
                default = getwd(),
                help    = paste("Path to save downloaded SRA files.",
                                "Default is CWD.")),
    make_option(c("-l", "--library_strategy"),
                type    = "character",
                default = NA,
                help    = paste("Filter for a specific library_strategy, for example:",
                                "ChIP-Seq, RNA-Seq, miRNA-Seq, ATAC-seq,",
                                "MeDIP-Seq, Bisulfite-Seq, WGS. Default is 'NA'."))
)
# Read options from command line
parser     <- OptionParser(usage = "%prog [options] SRA/SRP/SRX/SRR",
                           option_list = opt.list)
arguments  <- parse_args(parser, positional_arguments = c(1, Inf))
opt        <- arguments$options
accessions <- arguments$args

suppressPackageStartupMessages(library("SRAdb"))

# Check if SRA databases exists and download if needed
message("Checking if database exists...")
if (!file.exists(opt$database)) {
    message("Database file not found.")
    db.path  <- system(paste("dirname", opt$database), intern = TRUE)
    db.file  <- system(paste("basename", opt$database), intern = TRUE)
    if (!dir.exists(db.path)) {
        dir.create(db.path)
    }
    message(paste("Downloading database to", file.path(db.path, db.file)))
    sql.file <- getSRAdbFile(destdir  = db.path,
                             destfile = paste0(db.file, ".gz"))
} else {
    message(paste("Database file found at", opt$database))
    sql.file <- opt$database
}

### create a connection for later queries
sra_con <- dbConnect(SQLite(), sql.file)


suppressPackageStartupMessages(library("GEOmetadb"))
if(!file.exists(opt$geometadb)) getSQLiteFile()
### Check if GEO databases exists and download if needed
message("Checking if GEOmetadb database exists...")
if (!file.exists(opt$geometadb)) {
    message("Database file not found.")
    db.path  <- system(paste("dirname", opt$geometadb), intern = TRUE)
    db.file  <- system(paste("basename", opt$geometadb), intern = TRUE)
    if (!dir.exists(db.path)) {
        dir.create(db.path)
    }
    message(paste("Downloading database to", file.path(db.path, db.file)))
    sql.file <- getSQLiteFile(destdir  = db.path,
                              destfile = paste0(db.file, ".gz"))
} else {
    message(paste("Database file found at", opt$geometadb))
    sql.file <- opt$geometadb
}
geo_con <- dbConnect(SQLite(), opt$geometadb)

### Create output directory if it doesn't exist
message("Checking if output directory exists...")
if (!dir.exists(opt$outdir)) {
    message(paste("Creating output directory at", opt$outdir))
    dir.create(opt$outdir)
}

myGetAddress <- function(accession, sra_con,
                         fileType = 'sra', srcType='ftp') {
  return(listSRAfile(in_acc=accession, sra_con = sra_con,
         fileType = fileType, srcType=srcType)[,5])
}

myGetSampleAlias <- function(srs, sra_con) {
  return (unlist(dbGetQuery(sra_con, paste("SELECT sample_alias, scientific_name ",
    "FROM sample where sample_accession='", srs, "' ", sep="") ) ) )
}

myGetExperimentAlias <- function(srx, sra_con) {
  return (unlist(dbGetQuery(sra_con, paste("SELECT experiment_alias ",
    "FROM experiment WHERE experiment_accession='", srx, "'", sep="") ) ) )
}

myGetGSMTilte <- function(gsm, geo_con){
  return (unlist(dbGetQuery(geo_con,
    paste("SELECT title FROM gsm WHERE gsm='", gsm, "'", sep="") ) ) )
}

isGSM <- function(s) { return(substr(s, 1,3)=="GSM")}

for (accession in accessions) {
  conversion <- sraConvert(in_acc= accession, sra_con = sra_con,
      out_type = c("sra", "submission", "study", "sample", "experiment", "run") )

  conversion$ftp <- apply(matrix(conversion$run, ncol=1), MAR=1,
      FUN = myGetAddress, sra_con = sra_con, fileType = 'sra', srcType='ftp')

  conversion$experiment_alias <- apply(matrix(conversion$experiment, ncol=1),
      MAR=1, FUN=myGetExperimentAlias, sra_con=sra_con)

  conversion <- cbind(conversion, t(apply(matrix(conversion$sample, ncol=1),
      MAR=1, FUN=myGetSampleAlias, sra_con=sra_con)))

  info <- as.data.frame(matrix(NA, nrow(conversion), 5))
  names(info) <- c("library_name", "library_strategy", "library_layout",
      "sample_attribute", "title")

  for (i in 1:nrow(conversion)) {
    run_accession <- conversion$run[i]
    info[i, 1:4] <- dbGetQuery(sra_con, paste("SELECT library_name, library_strategy,
      library_layout, sample_attribute FROM sra WHERE run_accession='",
      run_accession, "'", sep=""))
    if (isGSM(conversion$sample_alias[i])) {
      info[i, 5] <- myGetGSMTilte(conversion$sample_alias[i], geo_con)
    }
  }

  conversion <- cbind(conversion, info)

  ### filter by library_strategy if specified
  if (!is.na(opt$library_strategy))
    conversion <- conversion[conversion$library_strategy == opt$library_strategy, ]

  if (nrow(conversion) == 0) {
    message(paste("No records for accession", accession,
                  "and library_strategy", opt$library_strategy, "found."))
  } else {
    outfile <- paste(opt$outdir, "/", accession, ".csv", sep="")
    message(paste("Writing records to:", outfile))
    write.csv(conversion, file=outfile)
  }
}
