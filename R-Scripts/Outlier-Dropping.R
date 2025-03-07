####
#### Outlier-Dropping.R: Determine which predictors' outliers are most
#### detrimental to auto-coder performance (see https://is.gd/ClassifieR#step-5)
####
#### For help, open a shell (terminal) in this directory and run:
####    Rscript Outlier-Dropping.R --help
####


##Check working directory
if (!endsWith(getwd(), "R-Scripts")) {
  stop("Change working directory to R-Scripts before running Outlier-Dropping.R")
}

# Parse command-line arguments ------------------------------------------------
suppressPackageStartupMessages(library(optparse))
source("Rscript-Opts.R")
optionList <- c(commonOpts,
                unlist(noncommonOpts$OutlierDropping),
                unlist(noncommonOpts$RunUMSOutlierDropping))
helptext <- paste0("   ", c(
  "Determine which predictors' outliers are most detrimental to auto-coder",
  "performance.\n\n",
  "See https://is.gd/ClassifieR#step-5."))
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
seed2 <- opt$seed2
wilcoxAlpha <- opt$`wilcox-alpha`
setupOnly <- !is.null(opt$`debug-setup`)
maxAttempts <- opt$`max-attempts`
fstListSize <- opt$`fst-list-size`
bestDrop <- opt$`best-drop-cols`
logCSV <- opt$`drop-log-csv`


# Subroutines -----------------------------------------------------------

##Find the outliers that most hurt classifier performance by calculating the
##  correlation between resamples' test performance and proportion of outliers
outlierCors <- function(fstList, dropped=NULL) {
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(purrr))

  ##Get resample info for a single forest: a dataframe with a row for each
  ##  resample, with columns for each resample's AccAUC and proportion of
  ##  outliers for each measure (excluding already-dropped outliers)
  resampInfo <- function(fst, dropped) {
    ##Get outlier data for each token
    fst$trainingData %>%
      select(ends_with("_Outlier"), -!!dropped) %>%
      ##Add resample identity for each token
      mutate(rowIndex = 1:n(), .before=1) %>%
      left_join(fst$pred %>%
                  select(rowIndex, Resample),
                by="rowIndex") %>%
      ##Remove tokens not selected in any resample (possible w/ bootstrap)
      filter(!is.na(Resample)) %>%
      ##Get proportions of outliers for each resample
      group_by(Resample) %>%
      summarise(across(ends_with("_Outlier"), mean), .groups="drop") %>%
      ##Add performance results for each resample
      left_join(fst$resample %>%
                  select(Resample, AccAUC),
                by="Resample")
  }

  ##Get the correlation of each measure to AccAUC, as a tibble with the largest
  ##  negative outlier first
  ret <-
    fstList %>%
    map_dfr(resampInfo, dropped=dropped) %>%
    ##Get correlation matrix
    select(AccAUC, ends_with("_Outlier")) %>%
    cor()
  ##Reduce correlation matrix to a vector, keeping just AccAUC row, and deleting
  ##  AccAUC's correlation with itself
  ret <- ret[1, -1] %>%
    ##Get tibble sorted by correlation
    tibble(Measure = names(.),
           TestCor = .) %>%
    arrange(TestCor)

  ret
}

##Run classifier with outlier columns included
outlierCls <- function(formula, data, UMS=strategy, smry=clsSummaryFunc,
                       tryDrop, dropped, attempt,
                       seedListSeed=seed2, clsSeed=seed,
                       tuneParams=clsTuneGrid, nTrees=numTrees,
                       stopAfterSetup=setupOnly,
                       savePrefix=outputPrefix, saveDir=outputDir) {
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(tidyr))

  ##Set up variables
  ##If tryDrop is null, then it's just the initial run to get the first ranking of outliers
  initRun <- is.null(tryDrop)
  ##numDropped is 0 for non-initial runs, or trying to drop the first outlier
  numDropped <- length(dropped)

  ##Prevent accidentally re-trying a measure that's already been dropped (only
  ##  applicable if outliers have already been dropped)
  if (!initRun && numDropped > 0 && tryDrop %in% dropped) {
    stop("Measure ", tryDrop, " has already been dropped.")
  }
  ##Attempt must be 1, 2, or 3 (unless tryDrop is null, when attempt can be anything)
  if (!initRun) {
    stopifnot(attempt %in% 1:3)
  }

  ##Print header
  if (initRun) {
    dropRound <- "0"
    printHeader("Outliers drop 0: Initial run to establish outlier ranking")
  } else {
    dropRound <- paste(numDropped+1, attempt, sep=".")
    printHeader(paste0("Outliers drop ", dropRound, ": ", tryDrop))
    if (numDropped > 0) {
      cat("Already been dropped:", paste(dropped, collapse=", "), fill=TRUE)
    }
  }

  ##Drop tokens with outliers (both established and new)
  if (initRun) {
    data <- umsData(data, UMS, outlierCols=TRUE)
  } else {
    data <- umsData(data, UMS, outlierCols=TRUE, dropCols=c(tryDrop, dropped))
  }

  ##Display data
  printData(data, UMS, outlierCols=TRUE)

  ##Optionally skip running classifier list (for debugging)
  if (stopAfterSetup) {
    cat("Skipping classifier run", fill=TRUE)
    return(NULL)
  }

  ##Run classifier list
  suppressPackageStartupMessages(library(caret))
  suppressPackageStartupMessages(library(foreach))
  suppressPackageStartupMessages(library(doParallel))

  ##Register parallel cluster
  registerDoSEQ()

  ##Set parallel seeds
  set.seed(seedListSeed)
  ##For bootstrap resampling with k resamples, trainControl(seeds) must be a list
  ##  of k+1 numeric vectors, each containing M integers (M = number of models being
  ##  evaluated within each call to train())
  seedList <- replicate(fstListSize,
                        replicate(26, sample.int(10000, 1), simplify=FALSE),
                        simplify=FALSE)
  cl <- makeCluster(fstListSize)
  registerDoParallel(cl)
  clusterSetRNGStream(cl, clsSeed)

  ##Train forest list
  cat("Running classifier list...", fill=TRUE)
  cls <- foreach(seed = seedList,
                 .errorhandling="pass",
                 .packages=c("caret","dplyr")) %dopar% {
                   train(formula,
                         data=data,
                         method="ranger",
                         ##caret/ranger default settings
                         importance="impurity",
                         num.trees=nTrees,
                         tuneGrid=tuneParams,
                         trControl=trainControl(method="boot",
                                                number=25,
                                                sampling=NULL,
                                                seeds=seed,
                                                savePredictions=TRUE,
                                                summaryFunction=smry,
                                                returnResamp="all",
                                                classProbs=TRUE))
                 }

  ##Stop the parallel cluster
  stopCluster(cl)

  ##Add name attribute
  clsName <- paste0(savePrefix, "_UMS", UMS, "_round", dropRound, ".Rds")
  attr(cls, "name") <- clsName

  ##Optionally save classifier list as Rds
  cat("Classifiers finished running", fill=TRUE)
  if (!is.null(saveDir)) {
    cat("Saving classifier list as", clsName, fill=TRUE)
    saveRDS(cls, paste0(saveDir, "/", clsName))
  } else {
    cat("Not saving classifier list", fill=TRUE)
  }

  ##Return classifier
  cls
}

##Get resampled performance information (flattened if it's a classifier list)
clsResamples <- function(x) {
  suppressPackageStartupMessages(library(purrr))
  suppressPackageStartupMessages(library(dplyr))

  ##Check that x is either a classifier or a classifier list
  is_train <- !is.null(x) && "train" %in% class(x)
  is_train_list <- map_lgl(x, ~ "train" %in% class(.x)) %>% all()
  if (!is_train && !is_train_list) {
    stop("x must be a classifier or classifier list run using caret::train()")
  }

  if (is_train_list) {
    resamp <- x %>%
      map_dfr("resample", .id="Model")
  } else {
    resamp <- x$resample
  }
  resamp %>%
    relocate(AccAUC, Accuracy, AUC, .before=1)
}


# Initial run -----------------------------------------------------------------
suppressPackageStartupMessages(library(dplyr))

##Display headers & define classifier info
source("UMS-Utils.R")
printHeader("Outlier dropping", "=")
printHeaderUMS(strategy)
##Unlike other scripts, keep outlier columns in clsData and exclude them via
##  model formula (i.e., Rpresent ~ . - race - Outlier1 - Outlier2 - ...).
##This keeps the full dataframe (outlier columns included) in the output
##  classifier object, which in turn makes it easier to assess the impact of
##  outlier columns on classifier performance (using outlierCors())
##umsData() is used to apply UMS within outlierCls()
clsData <- readRDS(inputPath)
clsData <- clsData[clsData$HowCoded=="Hand",]
##Get model formula that excludes all outliers in the UMS'd version of the data
##  (i.e., Rpresent ~ . - race - Outlier1 - Outlier2 - ...). Need to use
##  umsData() here because of column-dropping UMSs (i.e., if the UMS drops Col1,
##  clsForm should *not* have `- Col1_Outlier`).
clsForm <- umsFormula(umsData(clsData, strategy, outlierCols=TRUE),
                      strategy, outlierCols=TRUE)
clsTuneGrid <- umsTuneGrid(strategy, tuneGrid)
clsSummaryFunc <- umsSummaryFunc(strategy)

##Get output prefix
if (is.null(outputPrefix)) {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  fileArg <- "--file="
  scriptName <- sub(fileArg, "", cmdArgs[grep(fileArg, cmdArgs)])
  outputPrefix <- sub(".R", "", scriptName)
}

##Initial run
drop0 <- outlierCls(clsForm, clsData, tryDrop=NULL, dropped=character(0L))
if (setupOnly) {
  stop("Setup complete")
}

##Log results
drop0perf <- drop0 %>%
  clsResamples() %>%
  summarise(across(AccAUC, list(mean=mean, sd=sd)))
dropLog <- cbind(
  data.frame(round = 0, dropped = NA, tryDrop = NA,
             numDropped = 0, attempt = NA, classifierName = attr(drop0, "name")),
  drop0perf,
  data.frame(wilcoxW = NA, wilcoxP = NA, dropSuccess = NA)
)
if (!is.null(logCSV)) {
  write.csv(dropLog, logCSV, row.names=FALSE)
}




# Main loop ---------------------------------------------------------------
##Set looping variables to initial values
baseline <- drop0
tryDrop <- NULL
dropped <- character(0L)
attempt <- 1

##Loop over procedure to attempt to add a new outlier
##Stop if we've maxed out attempts (or, if maxAttempts is 0, skip the loop entirely)
while (attempt <= maxAttempts) {
  ##Determine which outliers most negatively impact performance (AccAUC) in the
  ##  baseline classifier list
  negOutliers <- outlierCors(baseline, dropped=dropped)
  ##Try dropping the measure corresponding to this attempt
  tryDrop <- negOutliers$Measure[attempt]

  ##Run a new classifier list
  new <- outlierCls(clsForm, clsData,
                    tryDrop=tryDrop, dropped=dropped, attempt=attempt)

  ##Store labels for models based on last outlier(s) dropped
  numDropped <- length(dropped)
  lastDropped <- ifelse(numDropped > 0,
                        dropped[numDropped],
                        "None")
  dropLabels <- gsub("_Outlier", "", c(lastDropped, tryDrop))

  ##Run Wilcoxon test to asses whether dropping tryDrop significantly improves performance
  ##Get performance (Accuracy times AUC) results for each classifier within each Dropped level
  perf <-
    ##Put resamples for both baseline and new into a dataframe
    list(baseline, new) %>%
    set_names(dropLabels) %>%
    map_dfr(clsResamples, .id="Dropped") %>%
    ##Calculate each model's performance
    group_by(Dropped, Model) %>%
    summarise(across(AccAUC, mean, na.rm=TRUE), .groups="drop") %>%
    ##Order Dropped chronologically (since Wilcoxon test is one-tailed)
    mutate(across(Dropped, factor, levels=dropLabels))
  ##Run Wilcoxon test
  outlierWilcox <- wilcox.test(AccAUC ~ Dropped, perf, alternative="less", exact=TRUE)
  ##Print test info
  wilcoxW <- outlierWilcox$statistic
  wilcoxP <- outlierWilcox$p.value
  dropSuccess <- wilcoxP < wilcoxAlpha
  cat("Wilcoxon rank sum exact test (Ha: mean < 0): W = ",
      wilcoxW, ", p-value ",
      if (!startsWith(format.pval(wilcoxP), "<")) "= ",
      format.pval(wilcoxP),
      sep="", fill=TRUE)

  ##Log info
  ##Format info for log
  if (numDropped==0) {
    logDropped <- "None"
  } else {
    logDropped <- gsub("_Outlier", "", paste(dropped, collapse=", "))
  }
  logTryDrop <- dropLabels[2]
  dropRound <- paste(numDropped+1, attempt, sep=".")
  ##Log
  newPerf <-
    perf %>%
    filter(Dropped==logTryDrop) %>%
    summarise(across(AccAUC, list(mean=mean, sd=sd)))
  dropLog <- rbind(
    dropLog,
    cbind(
      data.frame(round = dropRound, dropped = logDropped, tryDrop = logTryDrop,
                 numDropped = numDropped, attempt = attempt,
                 classifierName = attr(new, "name")),
      newPerf,
      data.frame(wilcoxW = wilcoxW, wilcoxP = wilcoxP, dropSuccess = dropSuccess))
  )

  ##Optionally write log to file
  if (!is.null(logCSV)) {
    write.csv(dropLog, logCSV, row.names=FALSE)
  }

  ##Check significance
  if (dropSuccess) {
    cat("Dropping this outlier significantly **improves** classifier performance", fill=TRUE)
    ##If significant, add tryDrop to list of dropped outliers, set new
    ##  classifier as baseline, and reset attempts counter
    dropped <- c(dropped, tryDrop)
    baseline <- new
    attempt <- 1
    ##Optionally write dropped outliers to file
    if (!is.null(bestDrop)) {
      writeLines(paste0(paste(dropped, collapse=","), "\n"), bestDrop)
    }
  } else {
    cat("Dropping this outlier **does not** significantly improve classifier performance", fill=TRUE)
    ##If nonsig, increment attempt
    attempt <- attempt + 1
  }

  ##Loop continues until attempt > maxAttempts
}


# After loop -----------------------------------------------------------------

##Report outcome
printHeader("All finished!")
if (length(dropped)==0) {
  cat("No outliers dropped", fill=TRUE)
} else {
  cat("The following outliers were dropped:", dropped, sep="\n", fill=TRUE)
}

##Save list of packages
pkg <- names(sessionInfo()$otherPkgs)
cat(pkg, file=pkglistPath, fill=TRUE, append=TRUE)
