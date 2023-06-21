####
#### Rscript-Opts.R: Utility script defining command-line arguments for
#### Run-UMS.R, Hyperparam-Tuning.R, Outlier-Dropping.R, and Session-Info.R.
#### 

suppressPackageStartupMessages(library(optparse))

##Options used by all 3 main scripts
commonOpts <- list(
  make_option(c("--ums", "-u"), type="character", default="0.0", 
              metavar="UMS-code", 
              help="Unfairness mitigation strategy (X.Y or X.Y.Z)"),
  make_option(c("--input-data", "-i"), type="character", metavar="path",
              default="../Input-Data/trainingData.Rds",
              help="Training data Rds file"),
  make_option(c("--seed", "-s"), type="integer", default=412412, metavar="seed",
              help="Random seed to set prior to running classifier"),
  make_option(c("--num-trees", "-n"), type="integer", default=500, 
              metavar="n",
              help="Number of trees per forest"),
  make_option(c("--output-dir", "-o"), type="character", metavar="path",
              default="../Outputs/Diagnostic-Files/Temp-Autocoders",
              help="Directory for saving classifier (NULL to skip saving)"),
  make_option(c("--output-prefix", "-p"), type="character", default="NULL",
              metavar="prefix",
              help="Prefix for saved classifier (NULL to use this script name)"),
  make_option(c("--pkglist-path", "-k"), type="character", metavar="path",
              default="../Outputs/Other/packages.tmp",
              help="tmpfile for saving package list")
)

##Options used by 1 or 2 main scripts
noncommonOpts <- list(
  ##Used exclusively by Run-UMS.R
  RunUMS = list(
    make_option(c("--drop-cols", "-d"), type="character", metavar="cols",
                default="NULL",
                help="Comma-separated list of columns for which rows with FALSE should be dropped (either string representation or file; NULL for none)")),
  
  ##Used exclusively by Hyperparam-Tuning.R
  HyperparamTuning = list(
    make_option(c("--model-status", "-m"), type="character", metavar="path",
                default="../Outputs/Diagnostic-Files/Model-Status/",
                help="Path to directory for saving tmpfiles indicating model status (NULL to skip saving)"),
    make_option(c("--best-tune-grid", "-b"), type="character", metavar="path",
                default="../Outputs/Other/Best-Params.csv",
                help="Path to csv for saving tuning parameters with best AccAUC (NULL to skip saving)")
  ),
  
  ##Used exclusively by Outlier-Dropping.R
  OutlierDropping = list(
    make_option(c("--seed2", "-S"), type="integer", default=302302,
                metavar="seed",
                help="Random seed to set prior to generating seed-list"),
    make_option(c("--wilcox-alpha", "-a"), type="double", default=.05,
                metavar="alpha-level",
                help="Alpha-level for Wilcoxon tests comparing performance based on\noutliers"),
    make_option(c("--debug-setup", "-g"), action="store_true",
                help="Stop before running classifiers (for debugging setup)"),
    make_option(c("--max-attempts", "-m"), type="integer", default=3,
                metavar="n",
                help="Number of consecutive unsuccessful outliers to stop after"),
    make_option(c("--fst-list-size", "-f"), type="integer", default=10,
                metavar="n",
                help="Number of forests per list"),
    make_option(c("--best-drop-cols", "-l"), type="character", metavar="path",
                default="../Outputs/Other/dropped-outliers.tmp",
                help="tmpfile for saving list of dropped outliers (NULL to skip saving)"),
    make_option(c("--drop-log-csv", "-c"), type="character", metavar="path",
                default="../Outputs/Other/Drop-Log.csv",
                help="csv for saving detailed log of dropped outliers (NULL to skip saving)")
  ),
  
  ##Used by Run-UMS.R and Outlier-Dropping.R
  RunUMSOutlierDropping = list(
    make_option(c("--tune-grid", "-t"), type="character", 
                metavar="dataframe-or-path",
                default="data.frame(mtry=13, splitrule=\"gini\", min.node.size=1)",
                help="Dataframe of hyperparameters")
  )
)

##Option used by Session-Info.R (similar to -k above but slightly different
##  help)
sinfoOpts <- list(
  make_option(c("--pkglist-path", "-k"), type="character", metavar="path",
              default="../Outputs/Other/packages.tmp",
              help="tmpfile where package list is stored")
)

##Display helptext if this script called directly
if ("--file=Rscript-Opts.R" %in% commandArgs()) {
  cat("Rscript-Opts.R:",
      paste0("   ", c(
        "Utility script defining command-line arguments for Run-UMS.R,",
        "Hyperparam-Tuning.R, Outlier-Dropping.R, and Session-Info.R.",
        "\n",
        "Not meant to be run directly with Rscript."
      )), fill=TRUE)
}
