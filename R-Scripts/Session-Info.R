####
#### Session-Info.R: Print session info based on packages listed in a textfile.
#### Useful for when a task is split across multiple independent R processes,
#### such as multiple Rscript calls in a shell script.
####
#### For help, open a shell (terminal) in this directory and run:
####    Rscript Session-Info.R --help
####

##Check working directory
if (!endsWith(getwd(), "R-Scripts")) {
  stop("Change working directory to R-Scripts before running Session-Info.R")
}

##Parse command line argument
suppressPackageStartupMessages(library(optparse))
source("Rscript-Opts.R")
optionList <- sinfoOpts
helptext <- paste0("   ", c(
  "Print session info based on packages listed in a textfile.",
  "",
  "Useful for when a task is split across multiple independent R processes,",
  "such as multiple Rscript calls in a shell script."))
opt_parser <- OptionParser(option_list=optionList, description=helptext)
opt <- parse_args(opt_parser)
filename <- opt$`pkglist-path`

##Read tmpfile
if (!file.exists(filename)) {
  stop("File ", filename, " does not exist.")
}
pkgList <- readLines(filename)

##Load packages
pkg <- unique(unlist(strsplit(pkgList, " ")))
libSilent <- function(x) {
  stopifnot(is.character(x))
  suppressPackageStartupMessages(library(x, character.only=TRUE))
}
invisible(lapply(pkg, libSilent))

##Print session info
source("UMS-Utils.R")
printHeader("Session info")
print(sessionInfo())
