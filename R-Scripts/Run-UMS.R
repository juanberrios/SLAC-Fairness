####
#### Run-UMS.R: Generate an auto-coder according to an unfairness mitigation
#### strategy
####
#### For help, open a shell (terminal) in this directory and run:
####    Rscript Run-UMS.R --help
####

##Check working directory
if (!endsWith(getwd(), "R-Scripts")) {
  stop("Change working directory to R-Scripts before running Run-UMS.R")
}

# Parse command line arguments --------------------------------------------
suppressPackageStartupMessages(library(optparse))
source("Rscript-Opts.R")
optionList <- c(commonOpts,
                unlist(noncommonOpts$RunUMS),
                unlist(noncommonOpts$RunUMSOutlierDropping))
helptext <- paste0("   ", c(
  "Generate an auto-coder according to an unfairness mitigation strategy"))
opt_parser <- OptionParser(option_list=optionList, description=helptext)
opt <- parse_args(opt_parser)

##Pull out arguments
##Common
strategy <- opt$ums
inputPath <- opt$`input-data`
seed <- opt$seed
numTrees <- opt$`num-trees`
outputDir <- opt$`output-dir`
outputDir <- if (is.null(outputDir) || outputDir=="NULL") NULL else outputDir
outputPrefix <- opt$`output-prefix`
outputPrefix <- if (is.null(outputPrefix) || outputPrefix=="NULL") NULL else outputPrefix
pkglistPath <- opt$`pkglist-path`
##Non-common
tuneGrid <- opt$`tune-grid`
dropCols <- opt$`drop-cols`

##Parse dropCols
if (is.null(dropCols) || dropCols=="NULL") {
  dropCols <- NULL
} else {
  if (file.exists(dropCols)) {
    dropCols <- readLines(dropCols)
  }
  dropCols <- strsplit(dropCols, ",")[[1]]
  dropCols <- dropCols[nchar(dropCols)>0]
}


# Main execution ----------------------------------------------------------

##Display headers & define classifier info
source("UMS-Utils.R")
printHeaderUMS(strategy)
clsData <- readRDS(inputPath)
clsData <- umsData(clsData, strategy, dropCols=dropCols)
printData(clsData, strategy)
clsForm <- umsFormula(clsData, strategy)
clsTuneGrid <- umsTuneGrid(strategy, tuneGrid)
clsSummaryFunc <- umsSummaryFunc(strategy)

##Run classifier
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(ranger))
cat("Running classifier...", fill=TRUE)
set.seed(seed)
##Separate for UMS 0.2
if (strategy=="0.2") {
  ##Get group levels
  gp <- clsData$race ##FYI: Modify if group isn't race
  if (!is.factor(gp)) {
    gp <- factor(gp)
  }
  levs <- setNames(levels(gp), levels(gp))

  suppressPackageStartupMessages(library(purrr))
  cls <- map(levs,
             ~ train(clsForm,
                     data=clsData[gp==.x,],
                     method="ranger",
                     importance="impurity",
                     num.trees=numTrees,
                     tuneGrid=clsTuneGrid,
                     trControl=trainControl(method="boot",
                                            number=25,
                                            sampling=NULL,
                                            savePredictions=TRUE,
                                            summaryFunction=clsSummaryFunc,
                                            returnResamp="all",
                                            classProbs=TRUE))
  )
##Single classifier for all other UMSs
} else {
  cls <- train(clsForm,
               data=clsData,
               method="ranger",
               importance="impurity",
               num.trees=numTrees,
               tuneGrid=clsTuneGrid,
               trControl=trainControl(method="boot",
                                      number=25,
                                      sampling=NULL,
                                      savePredictions=TRUE,
                                      summaryFunction=clsSummaryFunc,
                                      returnResamp="all",
                                      classProbs=TRUE))
}
cat("Classifier finished running", fill=TRUE)

##Save classifier as Rds
if (!is.null(outputDir)) {
  ##If no output prefix specified, get calling file name
  if (is.null(outputPrefix)) {
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    fileArg <- "--file="
    scriptName <- sub(fileArg, "", cmdArgs[grep(fileArg, cmdArgs)])
    outputPrefix <- sub(".R", "", scriptName)
  }
  ##Add UMS suffix and ".Rds"
  clsName <- paste0(outputPrefix, "_UMS", strategy, ".Rds")
  ##Save
  cat("Saving classifier as", clsName, fill=TRUE)
  saveRDS(cls, paste0(outputDir, "/", clsName))
} else {
  cat("Not saving classifier", fill=TRUE)
}

##Save list of packages
pkg <- names(sessionInfo()$otherPkgs)
cat(pkg, file=pkglistPath, fill=TRUE, append=TRUE)
