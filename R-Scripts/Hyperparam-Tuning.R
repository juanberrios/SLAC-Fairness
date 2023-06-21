####
#### Hyperparam-Tuning.R: Determine which hyperparameter settings optimize
#### auto-coder performance (see https://is.gd/ClassifieR#step-4)
####
#### For help, open a shell (terminal) in this directory and run:
####    Rscript Session-Info.R --help
####


##Parse command line arguments
suppressPackageStartupMessages(library(optparse))
source("Rscript-Opts.R")
optionList <- c(commonOpts, 
                unlist(noncommonOpts$HyperparamTuning))
helptext <- paste0("   ", c(
  "Determine which hyperparameter settings optimize auto-coder performance.",
  "",
  "See https://is.gd/ClassifieR#step-4."))
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
modStatusDir <- opt$`model-status`
paramsCSV <- opt$`best-tune-grid`

##Display headers & define classifier info
source("UMS-Utils.R")
printHeader("Hyperparameter tuning", "=")
printHeaderUMS(strategy)
clsData <- readRDS(inputPath)
clsData <- umsData(clsData[clsData$HowCoded=="Hand",], strategy)
printData(clsData, strategy)
clsForm <- umsFormula(clsData, strategy)
clsSummaryFunc <- umsSummaryFunc(strategy)

##Packages
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(ROCR))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))

##Register parallel cluster
cat("Registering cluster...", fill=TRUE)
##Define tuneList, whose elements are 1-row dataframes
tuneList <- expand.grid(mtry=seq(9, 17, by=2), 
                        splitrule=c("gini", "extratrees"), 
                        min.node.size=c(1, 5, 10), 
                        stringsAsFactors=FALSE) %>% 
  split(seq_len(nrow(.))) %>%
  set_names(sapply(., paste, collapse=" "))
##Reset any possible clusters
# numWorkers <- getDoParWorkers()
# cat("Number of parallel workers at script start:", numWorkers, "\n")
# if (numWorkers > 1) registerDoSEQ()
registerDoSEQ()
cl <- makeCluster(10)
registerDoParallel(cl)
##Set a parallel seed for reproducibility
clusterSetRNGStream(cl, seed)

##Define model status directory
if (!is.null(modStatusDir)) {
  statusPrefix <- paste0(modStatusDir, "/", strategy, "_TuneGrid")
}

cat("Running classifier list...", fill=TRUE)
##Train forest list
cls <-
  foreach(tuneVals = tuneList, 
          .packages=c("caret","dplyr")) %dopar% {
            ##Create tmpfile to signal that the model has started running
            file.create(paste(statusPrefix,
                              paste(tuneVals, collapse="&"), 
                              "begun.tmp", sep="_"))
            ##Call train()
            mod <- train(clsForm,
                         data=clsData,
                         method="ranger",
                         importance="impurity",
                         num.trees=numTrees,
                         tuneGrid=tuneVals,
                         trControl=trainControl(method="boot", 
                                                number=25,
                                                sampling=NULL,
                                                savePredictions=TRUE,
                                                summaryFunction=clsSummaryFunc,
                                                returnResamp="all",
                                                classProbs=TRUE))
            ##tmpfile: model has finished running
            file.create(paste(statusPrefix,
                              paste(tuneVals, collapse="&"), 
                              "completed.tmp", sep="_"))
            ##Return the classifier
            return(mod)
          } %>% 
  set_names(names(tuneList))


##Stop the parallel cluster
stopCluster(cl)
cat("Classifiers finished running", fill=TRUE)

##Display & optionally save optimal tuning params
##Get optimal tuning params
bestParams <- 
  cls %>% 
  ##Get resample element from each classifier object
  map_dfr("resample") %>% 
  ##Get classifier with greatest AccAUC
  group_by(mtry, splitrule, min.node.size) %>% 
  summarize(across(AccAUC, mean, na.rm=TRUE), .groups="drop") %>% 
  arrange(desc(AccAUC)) %>% 
  slice(1) %>%
  ##Just tuning grid
  select(-AccAUC)
##Display
disp <- 
  bestParams %>% 
  map_if(is.character, ~ paste0('"', .x, '"')) %>% 
  imap(~ paste(.y, "=", .x)) %>% 
  paste(collapse=", ") %>% 
  paste0("data.frame(", ., ")")
cat("Best tuning params:\n ", disp, fill=TRUE)
##Optionally save best tuning params
if (!is.null(paramsCSV)) {
  write.csv(bestParams, paramsCSV, row.names=FALSE)
  cat("Saving best params as", paramsCSV, fill=TRUE)
}

##If no output prefix specified, get calling file name
if (is.null(outputPrefix)) {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  fileArg <- "--file="
  scriptName <- sub(fileArg, "", cmdArgs[grep(fileArg, cmdArgs)])
  outputPrefix <- sub(".R", "", scriptName)
}
##Add name attribute
clsName <- paste0(outputPrefix, "_UMS", strategy, ".Rds")
attr(cls, "name") <- clsName

##Optionally save classifier list as Rds
if (!is.null(outputDir)) {
  cat("Saving classifier list as", clsName, fill=TRUE)
  saveRDS(cls, paste0(outputDir, "/", clsName))
} else {
  cat("Not saving classifier list", fill=TRUE)
}

##Save list of packages
pkg <- names(sessionInfo()$otherPkgs)
cat(pkg, file=pkglistPath, fill=TRUE, append=TRUE)
