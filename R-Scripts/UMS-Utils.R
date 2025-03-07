####
#### UMS-Utils.R: Utility functions for generating and analyzing auto-coders
#### that utilize unfairness mitigation strategies
####
#### For help, function descriptions below
####
#### NOTE: If using source(), either ensure working directory is the same as
#### this script's directory, or use source(keep.source=TRUE)
####


# Auto-coder generation -------------------------------------------------------

##These functions provide arguments to caret::train() according to different
##  unfairness mitigation strategies (UMSs)

##Some shared function arguments are:
##
##  Argument    Type        Description
##  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##  data        dataframe   Training data, with one row per token, one column
##                          for the dependent variable, one for the group
##                          characteristic, and multiple columns for predictors
##                          (aka features or acoustic measures)
##  strategy    character   UMS code (X.Y or X.Y.Z). Must be in
##                          ../Input-Data/UMS-List.txt
##  dependent   tidyselect* Column with variable to be auto-coded. Must be a
##                          two-level character vector or factor
##  group       tidyselect* Column with groups to assess fairness for. Must be a
##                          two-level character vector or factor
##  outlierCols logical     Account for outlier columns? If outlier-dropping,
##                          should be set to TRUE. More info under each function
##  umsTxt      character   Path from this script to file with two tab-separated
##                          columns: UMS and Description
## * See https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html

##Prepare data for auto-coding: apply unfairness mitigation strategy, only keep
##  necessary columns, and optionally drop outliers
umsData <- function(data,
                    strategy,
                    dependent=Rpresent,
                    group=race,
                    ##predictors (tidyselect): Set of columns that comprise
                    ##  predictor set (aka feature set)
                    predictors=TokenDur:absSlopeF0,
                    ##dropCols (character vector): Columns for which rows with
                    ##  FALSE should be dropped (i.e., for outlier-dropping).
                    ##  Rows are dropped before applying UMS
                    dropCols=NULL,
                    ##outlierCols (logical): Include outlier columns in output?
                    outlierCols=FALSE,
                    ##seed (numeric): Random seed to set before downsampling
                    seed=51233,
                    ##clsDir (character): Path to directory with precursor
                    ##  auto-coders, for UMSs 2.1.x
                    clsDir="../Outputs/Diagnostic-Files/Temp-Autocoders/",
                    ##varimpDir (character): Path to directory with precursor
                    ##  auto-coder variable importance data, for UMS 2.1.x
                    varimpDir="../Outputs/Other/",
                    ##normFile (character): Path to csv file with speaker
                    ##  normalization baselines, for UMS 3.1
                    normFile="../Input-Data/meanPitches.csv",
                    umsTxt="../Input-Data/UMS-List.txt") {
  # Implement strategies ------------------------------------------------------
  suppressWarnings(suppressMessages(library(dplyr)))
  ##Subroutine for core UMS functionality
  implementUMS <- function(x, UMS) {
    ##If not normalizing, drop Speaker column (if it exists)
    if (!startsWith(UMS, "3")) {
      x <- x %>%
        select(-any_of("Speaker"))
    }

    # Baseline/precursor UMSs -----------------------------------------------
    if (startsWith(UMS, "0")) {
      if (UMS %in% c("0.0", "0.1.1", "0.2")) {
        x <- x

      } else if (UMS %in% c("0.1.2")) {
        x <- x %>%
          select(-matches("(^F0|SlopeF0)"))

      } else if (UMS == "0.1.3") {
        x <- x %>%
          select({{dependent}}, {{group}}, matches("(^F0|SlopeF0)"))

      }

      ##Complete cases
      x <- filter(x, complete.cases(x))
      return(x)
    }


    # Downsampling UMSs -------------------------------------------------------
    if (startsWith(UMS, "1")) {
      ##Complete cases
      x <- filter(x, complete.cases(x))

      if (UMS == "1.1") {
        ##Get token count from smaller group
        nSmallerGroup <-
          x %>%
          count({{group}}) %>%
          filter(n==min(n)) %>%
          pull(n)

        ##Downsample by group to smaller amount
        set.seed(seed)
        x <- x %>%
          ##Downsample men's data so there's a race balance
          group_by({{group}}) %>%
          slice_sample(n=nSmallerGroup) %>%
          ungroup()

      } else if (UMS == "1.2") {
        ##Get token count from minority class
        nMinorityClass <-
          x %>%
          count({{dependent}}) %>%
          filter(n==min(n)) %>%
          pull(n)

        ##Downsample by class to minority amount
        set.seed(seed)
        x <- x %>%
          ##Downsample men's data so there's a race balance
          group_by({{dependent}}) %>%
          slice_sample(n=nMinorityClass) %>%
          ungroup()

      } else if (UMS %in% c("1.3.1", "1.3.2", "1.4", "1.5", "1.6")) {
        ##Establish class & group hierarchies (majority/minority class, larger/smaller group)
        gps <-
          x %>%
          count({{group}}, sort=TRUE) %>%
          mutate(ord = c("larger","smaller")) %>%
          pull({{group}}, name=ord) %>%
          as.list()
        classes <-
          x %>%
          count({{dependent}}, sort=TRUE) %>%
          mutate(ord = c("majority","minority")) %>%
          pull({{dependent}}, name=ord) %>%
          as.list()

        if (UMS == "1.3.1") {
          ##Get each group's percentage of majority class
          pctMajority <-
            x %>%
            group_by({{group}}) %>%
            summarise(MajPct = mean({{dependent}}==classes$majority)) %>%
            pull(MajPct, name={{group}})
          ratMajMinLarger <- pctMajority[gps$larger] / (1 - pctMajority[gps$larger])

          ##Determine which class to downsample for smaller group
          if (pctMajority[gps$smaller] > pctMajority[gps$larger]) {
            downsamp <- classes$majority
            ratio <- ratMajMinLarger
          } else {
            downsamp <- classes$minority
            ratio <- 1/ratMajMinLarger
          }

          ##Get token count for class-to-not-downsample for smaller group
          nSmallerReference <-
            x %>%
            filter({{group}}==gps$smaller,
                   {{dependent}}!=downsamp) %>%
            nrow()
          ##Get smaller group's data for class-to-downsample
          set.seed(seed)
          xSmallerDownsamp <-
            x %>%
            filter({{group}}==gps$smaller,
                   {{dependent}}==downsamp) %>%
            ##Downsample smaller group's class-to-downsample
            slice_sample(n=floor(ratio*nSmallerReference))

          ##Add larger group's x & rest of smaller group's data
          x <- x %>%
            filter(!({{group}}==gps$smaller & {{dependent}}==downsamp)) %>%
            rbind(xSmallerDownsamp)

        } else if (UMS == "1.3.2") {
          ##Get each group's percentage of majority class
          pctMajority <-
            x %>%
            group_by({{group}}) %>%
            summarise(MajPct = mean({{dependent}}==classes$majority)) %>%
            pull(MajPct, name={{group}})
          ratMajMinSmaller <- pctMajority[gps$smaller] / (1 - pctMajority[gps$smaller])

          ##Determine which class to downsample for larger group
          if (pctMajority[gps$smaller] > pctMajority[gps$larger]) {
            downsamp <- classes$minority
            ratio <- 1/ratMajMinSmaller
          } else {
            downsamp <- classes$majority
            ratio <- ratMajMinSmaller
          }

          ##Get token count for class-to-not-downsample for smaller group
          nLargerReference <-
            x %>%
            filter({{group}}==gps$larger,
                   {{dependent}}!=downsamp) %>%
            nrow()

          ##Get larger group's data for class-to-downsample
          set.seed(seed)
          xLargerDownsamp <-
            x %>%
            filter({{group}}==gps$larger,
                   {{dependent}}==downsamp) %>%
            ##Downsample larger group's class-to-downsample
            slice_sample(n=round(ratio*nLargerReference))

          ##Add smaller group's x & rest of larger group's data
          x <- x %>%
            filter(!({{group}}==gps$larger & {{dependent}}==downsamp)) %>%
            rbind(xLargerDownsamp)

        } else if (UMS == "1.4") {
          ##Get smaller group's class distribution
          smallerGp <-
            x %>%
            filter({{group}}==gps$smaller) %>%
            count({{dependent}}) %>%
            pull(n, name={{dependent}})

          ##Downsample data
          set.seed(seed)
          ##Select group-balanced subsamples by class
          xMin <-
            x %>%
            filter({{dependent}}==classes$minority) %>%
            group_by({{group}}) %>%
            slice_sample(n=smallerGp[classes$minority]) %>%
            ungroup()
          xMaj <-
            x %>%
            filter({{dependent}}==classes$majority) %>%
            group_by({{group}}) %>%
            slice_sample(n=smallerGp[classes$majority]) %>%
            ungroup()

          ##Combine subsamples
          x <- rbind(xMin, xMaj)

        } else if (UMS == "1.5") {
          ##Get minority class's group distribution
          minorityClass <-
            x %>%
            filter({{dependent}}==classes$minority) %>%
            count({{group}}) %>%
            pull(n, name={{group}})

          ##Downsample data
          set.seed(seed)
          ##Select class-balanced subsamples by group
          xSmaller <-
            x %>%
            filter({{group}}==gps$smaller) %>%
            group_by({{dependent}}) %>%
            slice_sample(n=minorityClass[gps$smaller]) %>%
            ungroup()
          xLarger <-
            x %>%
            filter({{group}}==gps$larger) %>%
            group_by({{dependent}}) %>%
            slice_sample(n=minorityClass[gps$larger]) %>%
            ungroup()

          ##Combine subsamples
          x <- rbind(xSmaller, xLarger)

        } else if (UMS == "1.6") {
          ##Get smallest class/group sample
          nSmallest <-
            x %>%
            count({{dependent}}, {{group}}) %>%
            pull(n) %>%
            min()

          ##Downsample data
          set.seed(seed)
          x <- x %>%
            group_by({{dependent}}, {{group}}) %>%
            slice_sample(n=nSmallest) %>%
            ungroup()

        }
      }

      return(x)
    }


    # Valid predictor selection UMSs ------------------------------------------
    if (startsWith(UMS, "2")) {
      if (startsWith(UMS, "2.1")) {
        ##Load purrr
        suppressWarnings(suppressMessages(library(purrr)))

        ##Get precursor file listings
        fileStem <- if_else(grepl("2\\.1\\.[1-3]", UMS), "0.1.1", "0.2")
        predFiles <- list(
          Cls = dir(clsDir, paste0("UMS", fileStem), full.names=TRUE),
          VarImp = dir(varimpDir, paste0("Var-Imp_UMS", fileStem), full.names=TRUE)
        )

        ##Handle bad number of files
        numFiles <- map(predFiles, length)
        ##Handle both files missing
        if (numFiles$Cls==0 && numFiles$VarImp==0) {
          stop("UMS ", UMS, " requires either a UMS ", fileStem,
               " classifier file in clsDir or a UMS ", fileStem,
               " variable importance file in varimpDir.")
        }
        ##Handle duplicate files
        if (numFiles$Cls > 1 || numFiles$VarImp > 1) {
          taMsg <- "There are"
          clsMsg <- if (numFiles$Cls > 1)
            paste(" multiple UMS", fileStem, "classifier files in clsDir")
          andMsg <- if (numFiles$Cls > 1 && numFiles$VarImp > 1)
            ", and"
          varImpMsg <- if (numFiles$VarImp > 1)
            paste(" multiple UMS", fileStem, "variable importance files in varimpDir")
          stop(taMsg, clsMsg, andMsg, varImpMsg)
        }

        ##Handle varImp
        if (numFiles$VarImp == 1) {
          ##Read
          varImpCsv <- read.csv(predFiles$VarImp)

          ##Check format
          if (fileStem == "0.1.1") {
            varImpCols <- c("Importance", "Measure")
          } else if (fileStem == "0.2") {
            gps <-
              x %>%
              pull({{group}}) %>%
              unique()
            varImpCols <- c(paste0("Importance_", gps), "Measure")

          }
          if (!identical(sort(colnames(varImpCsv)), sort(varImpCols))) {
            stop("UMS ", fileStem, " variable importance file needs to have columns ",
                 paste(varImpCols, collapse=","))
          }

          ##Define most important variables from varImpCsv if no varImpCls
          if (numFiles$Cls == 0) {
            impVars <- varImpCsv
          }
        }

        ##Handle varImpCls
        if (numFiles$Cls == 1) {
          ##Create varImpCls dataframe
          if (fileStem == "0.1.1") {
            varImpCls <-
              readRDS(predFiles$Cls) %>%
              pluck("finalModel", "variable.importance") %>%
              {tibble(Measure=names(.), Importance=.)}
          } else if (fileStem == "0.2") {
            suppressWarnings(suppressMessages(library(tidyr)))
            varImpCls <-
              readRDS(predFiles$Cls) %>%
              map_dfr(~ .x %>%
                        pluck("finalModel", "variable.importance") %>%
                        {tibble(Measure=names(.), value=.)},
                      .id="name") %>%
              pivot_wider(names_prefix="Importance_")
          }

          ##Warn if info isn't the same between varImpCls & varImpCsv
          if (numFiles$VarImp == 1 && !isTRUE(all.equal(as.data.frame(varImpCls), varImpCsv))) {
            warning("UMS ", fileStem, " variable importance file conflicts with ",
                    "variable importance in classifier file.\n  (Consider ",
                    "regenerating variable importance file.)\n  ",
                    "Using classifier file for variable importance.")
          }

          ##Define most important variables from varImpCls
          impVars <- varImpCls
        }

        ##Determine variables to toss
        if (fileStem == "0.1.1") {
          ##For UMS 2.1.[123], rank vars by importance (in predicting race)
          impVars <- impVars %>%
            mutate(Rank = rank(desc(Importance)))

          ##Get rank cutoff (top X%)
          rankCutoff <- floor(nrow(impVars) * case_when(UMS == "2.1.1" ~ 0.1,
                                                        UMS == "2.1.2" ~ 0.2,
                                                        UMS == "2.1.3" ~ 0.5))

          ##Select vars ranked rankCutoff or higher
          varsToToss <-
            impVars %>%
            filter(Rank <= rankCutoff) %>%
            pull(Measure)
        } else if (fileStem == "0.2") {
          ##For UMS 2.1.[45], calculate absolute difference in importance ranks
          impVars <- impVars %>%
            relocate(Measure, .before=1) %>%
            mutate(across(starts_with("Importance"), ~ rank(desc(.x)))) %>%
            mutate(RankDiff = abs(.[[2]] - .[[3]]))

          ##Get rank diff cutoff (p/N rank places)
          rankDiffCutoff <- floor(nrow(impVars) / case_when(UMS == "2.1.4" ~ 2,
                                                            UMS == "2.1.5" ~ 3))

          ##Select vars with rank diff at least rankDiffCutoff
          varsToToss <-
            impVars %>%
            filter(RankDiff >= rankDiffCutoff) %>%
            pull(Measure)
        }

        ##Toss out variables (including outlier columns if applicable)
        if (outlierCols) {
          varsToToss <- c(varsToToss, paste0(varsToToss, "_Outlier")) %>%
            intersect(colnames(x))
        }
        x <- x %>%
          select(-all_of(varsToToss))
        ##End if (startsWith(UMS, "2.1"))
      } else if (UMS == "2.2") {
        ##Without F0 measures
        x <- x %>%
          select(-matches("(^F0|SlopeF0)"))

      } else if (UMS == "2.3") {
        ##Without F0 point measures; the other F0 measures weren't very important
        ##  in UMS 0.1.1 (classifier of race)
        x <- x %>%
          select(-matches("^F0(min|max)"))

      }

      ##Complete cases
      x <- filter(x, complete.cases(x))
      return(x)
    }


    # Normalization UMS(s) ----------------------------------------------------
    if (startsWith(UMS, "3")) {
      ##Complete cases
      x <- filter(x, complete.cases(x))

      if (UMS == "3.1") {
        ##Check that Speaker is in data
        if (!("Speaker" %in% colnames(x))) {
          stop("UMS 3.1 requires Speaker column in x")
        }

        ##Check that precursor file exists & is in the correct format
        if (!file.exists(normFile)) {
          stop("UMS 3.1 requires a file with normalization baselines, but the file ", normFile, " does not exist.")
        }
        if (!endsWith(normFile, ".csv")) {
          stop("The normalization baseline file for UMS 3.1 must be in .csv format")
        }
        normBase <- read.csv(normFile)
        if (!all(c("Speaker", "MinPitch", "MaxPitch") %in% colnames(normBase))) {
          stop("The normalization baseline file for UMS 3.1 must have at least columns Speaker,MinPitch,MaxPitch")
        }

        ##Check that all speakers in x are in normBase
        xSpkrs <-
          x$Speaker %>%
          unique() %>%
          as.character()
        meanSpkrs <-
          normBase$Speaker %>%
          unique() %>%
          as.character()
        missingSpkrs <- setdiff(xSpkrs, meanSpkrs)
        if (length(missingSpkrs) > 0) {
          stop("The following speakers in data are missing from normalization baseline file:\n", paste(missingSpkrs, collapse=" "))
        }

        ##Normalize pitch according to baseline
        x <- x %>%
          left_join(normBase %>% select(Speaker, MinPitch, MaxPitch),
                    by="Speaker") %>%
          mutate(F0min = F0min - MinPitch,
                 F0max = F0max - MaxPitch)

      }

      ##Remove extra columns
      x <- x %>%
        select(-c(Speaker, MinPitch, MaxPitch))
      return(x)
    }
  }


  # Validate args -------------------------------------------------------------
  ##Validate strategy
  umsValidate(strategy, umsTxt=umsTxt)

  ##Validate columns
  stopIfNoCol(data, {{dependent}}, "dependent", "data")
  stopIfNoCol(data, {{group}}, "group", "data")
  stopIfNoCol(data, {{predictors}}, "predictors", "data")
  if (!is.null(dropCols)) {
    stopIfNoCol(data, {{dropCols}}, "dropCols", "data")
  }

  ##Check properties of data...
  ##  data is a data.frame
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  ##  dependent/group are binary
  nDep <- data %>% pull({{dependent}}) %>% n_distinct(na.rm=TRUE)
  if (nDep != 2) {
    stop("dependent must be binary; selected dependent has ", nDep, " levels.")
  }
  nGroup <- data %>% pull({{group}}) %>% n_distinct(na.rm=TRUE)
  if (nGroup != 2) {
    stop("group must be binary; selected group has ", nGroup, " levels.")
  }
  ##  there are Outlier columns (if outlierCols=TRUE)
  if (outlierCols && !any(endsWith(colnames(data), "_Outlier"))) {
    stop("outlierCols=TRUE, but there are no columns ending in _Outlier")
  }

  ##Validate strategy
  umsValidate(strategy, umsTxt=umsTxt)


  # Main execution ------------------------------------------------------------
  ##Optionally drop outliers
  if (!is.null(dropCols)) {
    data <- data %>%
      filter(if_all({{dropCols}}, `!`))
  }

  ##Drop unused columns (pre-UMS, to avoid drop_na() dropping too many rows)
  data <- data %>%
    select({{dependent}}, {{group}}, {{predictors}}, ends_with("_Outlier"),
           any_of("Speaker"))
  if (!outlierCols) {
    data <- data %>%
      select(-ends_with("_Outlier"))
  }

  ##Implement UMS
  if (!startsWith(strategy, "4")) {
    ##Simple UMS
    data <- implementUMS(data, strategy)
  } else {
    ##Combination UMSs
    ##Decode UMS code as 4.X.Y: X encodes predictor selection/normalization UMS,
    ##  Y encodes downsampling UMS
    components <- strsplit(strategy, ".", fixed=TRUE)[[1]]
    strat1 <- c("2.1.1", "2.2", "2.3", "3.1")[as.numeric(components[2])]
    strat2 <- c("1.3.1", "1.3.2")[as.numeric(components[3])]
    ##Implement
    data <- data %>%
      implementUMS(strat1) %>%
      implementUMS(strat2)
  }

  ##Return data
  data
}

##Specify model formula based on UMS (e.g., Rpresent ~ . - race)
umsFormula <- function(data,
                       strategy,
                       dependent=Rpresent,
                       group=race,
                       ##outlierCols (logical): Append " - [feature1]_Outlier
                       ##  - [feature2]_Outlier - ..." to formula?
                       outlierCols=FALSE,
                       umsTxt="../Input-Data/UMS-List.txt") {
  suppressWarnings(suppressMessages(library(dplyr)))

  ##Validate columns
  stopIfNoCol(data, {{dependent}}, "dependent", "data")
  stopIfNoCol(data, {{group}}, "group", "data")

  ##Validate strategy
  umsValidate(strategy, umsTxt=umsTxt)

  ##Get base formula
  dep <- deparse(substitute(dependent))
  gp <- deparse(substitute(group))
  if (strategy %in% c("0.1.1", "0.1.2", "0.1.3")) {
    ##e.g., race ~ . - Rpresent
    formula <- paste(gp, "~ . -", dep)
  } else {
    ##e.g., Rpresent ~ . - race
    formula <- paste(dep, "~ . -", gp)
  }

  ##If outliers, add " - Outlier1 - Outlier2 - ..."
  ##This keeps outlier columns in classifier's trainingData element, but
  ##  excludes outliers from predictor set
  if (outlierCols) {
    outCols <-
      data %>%
      select(ends_with("_Outlier")) %>%
      colnames()
    if (length(outCols)==0) {
      warning("Data does not have any columns ending in _Outlier")
    } else {
      formula <- paste(formula, "-", paste(outCols, collapse=" - "))
    }
  }

  ##Convert to formula
  as.formula(formula)
}

##Handle tuning parameter grid input and override mtry for UMS 0.1.3
umsTuneGrid <- function(strategy,
                        ##tuneGrid: 1-row dataframe, as either a dataframe
                        ##  object, a path to a .csv, or a string representation
                        ##  like "data.frame(..."
                        tuneGrid,
                        umsTxt="../Input-Data/UMS-List.txt") {
  ##Parse tuneGrid
  if (!is.data.frame(tuneGrid) && !is.character(tuneGrid)) {
    stop("tune-grid must be a dataframe or character")
  }
  ##Handle character tuneGrid (can be a path or a string rep of a df)
  if (is.character(tuneGrid)) {
    if (length(tuneGrid) != 1) {
      stop("tune-grid must have length 1")
    }
    ##Attempt to read from path
    if (file.exists(tuneGrid)) {
      tryCatch(tuneGrid <- read.csv(tuneGrid),
               error = function(e) {
                 if (conditionMessage(e) == "no lines available in input") {
                   stop("File ", tuneGrid, " is empty")
                 } else {
                   stop("Unknown error in read.csv(\"",tuneGrid, "\")")
                 }
               }
      )
      ##If no path, attempt to parse
    } else {
      suppressWarnings(suppressMessages(library(dplyr))) ##In case tuneGrid uses tibble()
      tryCatch(tuneGrid <- eval(str2expression(tuneGrid)),
               error = function(e) {
                 stop("tune-grid incorrectly specified (see below)\n",
                      conditionMessage(e))
               })
    }
    ##Enforce formatting
    if (!is.data.frame(tuneGrid) || nrow(tuneGrid) != 1) {
      stop("tune-grid must be a dataframe with 1 row")
    }
    if (!identical(sort(colnames(tuneGrid)), c("min.node.size","mtry","splitrule"))) {
      stop("tune-grid must have 3 columns: min.node.size, mtry, splitrule")
    }
    if (!is.numeric(tuneGrid$mtry)) {
      stop("mtry must be numeric")
    } else if (tuneGrid$mtry < 1) {
      stop("mtry must be positive")
    }
    if (!is.character(tuneGrid$splitrule) || !(tuneGrid$splitrule %in% c("gini","extratrees"))) {
      stop("splitrule must be either \"gini\" or \"extratrees\", not ", tuneGrid$splitrule)
    }
    if (!is.numeric(tuneGrid$min.node.size)) {
      stop("min.node.size must be numeric")
    } else if (tuneGrid$min.node.size < 1) {
      stop("min.node.size must be positive")
    }
  }

  ##Override mtry for 0.1.3 (just has 4 predictors)
  umsValidate(strategy, umsTxt=umsTxt)
  if (strategy =="0.1.3") {
    tuneGrid$mtry <- 2
  }

  ##Ensure tuneGrid isn't a tibble (will cause caret::train() to error out)
  if ("tbl_df" %in% class(tuneGrid)) {
    tuneGrid <- as.data.frame(tuneGrid)
  }

  tuneGrid
}

##Specify classifier summary function based on UMS
umsSummaryFunc <- function(strategy,
                           umsTxt="../Input-Data/UMS-List.txt") {
  ##Validate strategy
  umsValidate(strategy, umsTxt=umsTxt)

  if (strategy %in% c("0.1.1", "0.1.2", "0.1.3")) {
    summaryrace
  } else {
    summaryRpresent
  }
}


# Classifier summary functions --------------------------------------------
##Not meant to be used directly---formatted specifically for the summaryFunction
##  arg of caret::trainControl()
##See https://topepo.github.io/caret/model-training-and-tuning.html#metrics

##Dependent = Rpresent
summaryRpresent <- function(data, lev=NULL, model=NULL, obsCol="obs",
                            returnDF=FALSE) {
  suppressWarnings(suppressMessages(library(ROCR)))
  suppressWarnings(suppressMessages(library(dplyr)))
  suppressWarnings(suppressMessages(library(tidyr)))
  suppressWarnings(suppressMessages(library(purrr)))
  suppressWarnings(suppressMessages(library(magrittr)))

  if (!is.null(obsCol) && obsCol!="obs") {
    data <- data %>% rename("obs" = obsCol)
  }

  preds <- with(data, prediction(Present, obs))
  auc <-
    preds %>%
    performance(measure="auc") %>%
    attr("y.values") %>%
    extract2(1)

  cutoff <-
    preds %>%
    performance(measure="acc") %>%
    attributes() %>%
    magrittr::extract(c("x.values", "y.values")) %>%
    map_dfc(extract2, 1) %>%
    filter(y.values==max(y.values)) %>%
    pull(x.values) %>%
    median()

  data <- data %>%
    mutate(PredCutoff = if_else(Present >= cutoff, "Present", "Absent"),
           PredHalf = if_else(Present >= 0.5, "Present", "Absent"),
           AccCutoff = PredCutoff==obs,
           AccHalf = PredHalf==obs)

  class_means <-
    data %>%
    group_by(obs) %>%
    summarise(across(starts_with("Acc"), mean)) %>%
    pivot_longer(starts_with("Acc"), names_to="Type", values_to="Acc") %>%
    mutate(Type = paste0(Type, "_", obs)) %>%
    select(Type, Acc) %>%
    pivot_wider(names_from=Type, values_from=Acc) %>%
    rename_with(~ gsub("Acc", "ClassAccuracy", .x)) %>%
    rename_with(~ gsub("Cutoff", "", .x))

  ret <- data.frame(AccAUC = mean(data$AccCutoff)*auc,
                    AUC = auc,
                    Accuracy = mean(data$AccCutoff),
                    BestCutoff = cutoff,
                    AccuracyHalf = mean(data$AccHalf)) %>%
    cbind(class_means)

  if (!returnDF) ret <- ret %>% map_dbl(I)
  return(ret)
}

##Dependent = race
summaryRender <- function(data, lev=NULL, model=NULL, obsCol="obs",
                          returnDF=FALSE) {
  suppressWarnings(suppressMessages(library(ROCR)))
  suppressWarnings(suppressMessages(library(dplyr)))
  suppressWarnings(suppressMessages(library(tidyr)))
  suppressWarnings(suppressMessages(library(purrr)))
  suppressWarnings(suppressMessages(library(magrittr)))

  if (!is.null(obsCol) && obsCol!="obs") {
    data <- data %>% rename("obs" = obsCol)
  }

  preds <- with(data, prediction(white, obs))
  auc <-
    preds %>%
    performance(measure="auc") %>%
    attr("y.values") %>%
    extract2(1)

  cutoff <-
    preds %>%
    performance(measure="acc") %>%
    attributes() %>%
    magrittr::extract(c("x.values", "y.values")) %>%
    map_dfc(extract2, 1) %>%
    filter(y.values==max(y.values)) %>%
    pull(x.values) %>%
    median()

  data <- data %>%
    mutate(PredCutoff = if_else(white >= cutoff, "white", "black"),
           PredHalf = if_else(white >= 0.5, "white", "black"),
           AccCutoff = PredCutoff==obs,
           AccHalf = PredHalf==obs)

  class_means <-
    data %>%
    group_by(obs) %>%
    summarise(across(starts_with("Acc"), mean)) %>%
    pivot_longer(starts_with("Acc"), names_to="Type", values_to="Acc") %>%
    mutate(Type = paste0(Type, "_", obs)) %>%
    select(Type, Acc) %>%
    pivot_wider(names_from=Type, values_from=Acc) %>%
    rename_with(~ gsub("Acc", "ClassAccuracy", .x)) %>%
    rename_with(~ gsub("Cutoff", "", .x))

  ret <- data.frame(AccAUC = mean(data$AccCutoff)*auc,
                    AUC = auc,
                    Accuracy = mean(data$AccCutoff),
                    BestCutoff = cutoff,
                    AccuracyHalf = mean(data$AccHalf)) %>%
    cbind(class_means)

  if (!returnDF) ret <- ret %>% map_dbl(I)
  return(ret)
}


# Fairness measurement --------------------------------------------------------
##These functions provide measurements to quantify an auto-coder's fairness

##Some shared function arguments are:
##
##  Argument    Type        Description
##  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##  x           object      Either an auto-coder to measure fairness for, or
##                          a predictions dataframe created by
##                          cls_fairness(x, output='pred')
##  group       tidyselect* Column with groups to assess fairness for. Must be a
##                          two-level character vector or factor
##  cutoffBy    character   Column of predictions dataframe defining groups of
##                          observations within which to find optimal
##                          classification thresholds. Default ("Resample")
##                          adjusts thresholds within each resample
##  unResample  logical     If there are multiple predictions per token in the,
##                          training data, TRUE (default) = adjust reported
##                          token counts (for chi-square or confusion matrix) to
##                          reflect the original training data size; FALSE =
##                          report total number of predictions
## * See https://dplyr.tidyverse.org/articles/programming.html#tidy-selection

##Investigate classifier fairness (meant for interactive use)
cls_fairness <- function(x,
                         ##output (character): Desired output.
                         ##  "dataframe": accuracy by group
                         ##  "chisq": homogeneity test for differences in
                         ##      accuracy between groups
                         ##  "cm": confusion matrices for each group
                         ##  "pred": predictions dataframe
                         output=c("dataframe","chisq","cm","pred"),
                         ##byClass (logical): For "dataframe" or "cm", report
                         ##  class accuracies (TRUE, default) or overall
                         ##  accuracy (FALSE)?
                         byClass=FALSE,
                         group=race,
                         cutoffBy=c("Resample","race"),
                         unResample=TRUE) {
  suppressWarnings(suppressMessages(library(dplyr)))
  suppressWarnings(suppressMessages(library(rlang)))

  # Subroutine: Get predictions dataframe -------------------------------------
  getPred <- function(cls) {
    # Check arguments ---------------------------------------------------------
    ##Check that x$pred exists
    rawPred <- cls$pred
    if (is.null(rawPred)) {
      stop("x is missing predictions; re-run caret::train() with ",
           "savePredictions=TRUE")
    }

    ##Check group column
    trainData <- ungroup(x$trainingData)
    stopIfNoCol(trainData, {{group}}, "group", "classifier's training data")
    groups <- unique(pull(trainData, {{group}}))
    numGroups <- length(groups)
    if (numGroups != 2) {
      groupsFormatted <- ifelse(numGroups > 6, c(groups[1:6], "..."), groups) %>%
        paste(collapse=" ")
      stop("group column (", as_name(enquo(group)), ") has ", numGroups,
           " groups (", groupsFormatted, ").",
           "Only two-group situations are supported")
    }

    ##Check binary classifier
    classes <- sort(unique(trainData$.outcome))
    numClasses <- length(classes)
    if (numClasses != 2) {
      classesFormatted <- ifelse(numClasses > 6, c(classes[1:6], "..."), classes) %>%
        paste(collapse=" ")
      stop("x has ", numClasses, " classes (", classesFormatted, ").",
           "Only two-class classifiers are supported")
    }

    # Calculate predictions ---------------------------------------------------
    ##Get classes
    nonProbClass <- classes[1]
    probClass <- classes[2]

    ##Add group from trainData to pred
    pred <-
      rawPred %>%
      ##Rename columns
      rename(Actual = obs, ProbOfProbClass = !!probClass) %>%
      ##Add group col from trainData
      left_join(trainData %>%
                  mutate(rowIndex = row_number()) %>%
                  ##We just need rowIndex & group
                  select(rowIndex, {{group}}),
                by="rowIndex")

    ##Get predictions, all with the same cutoff
    if (is.null(cutoffBy) || all(cutoffBy == "")) {
      ##Don't adjust cutoffs by group/resample
      pred <- pred %>%
        select(rowIndex, Actual, Predicted = pred,
               ##Rename ProbOfProbClass
               "Prob{probClass}" := ProbOfProbClass, {{group}})
      ##Add pred class attribute
      class(pred) <- c("pred", class(pred))

      return(pred)
    }

    ##Optionally adjust predictions for group/resample-specific cutoffs
    if (!is.null(cutoffBy) && all(cutoffBy != "")) {
      ##Check that conditions allow for specific cutoffs
      if (!(all(cutoffBy %in% colnames(pred)))) {
        missingCutoff <- paste(setdiff(cutoffBy, colnames(pred)), collapse=",")
        stop("Cannot calculate cutoffs by ", missingCutoff, "; missing from ",
             "x$pred")
      }

      suppressWarnings(suppressMessages(library(ROCR)))
      suppressWarnings(suppressMessages(library(purrr)))

      ##Calculate best cutoffs for each group
      cutoffs <-
        pred %>%
        group_by(pick(!!cutoffBy)) %>%
        group_modify(~ with(.x, prediction(ProbOfProbClass, Actual)) %>%
                       performance(measure="acc") %>%
                       ##Pull out Accuracy and Cutoff as a dataframe
                       attributes() %>%
                       `[`(c("x.values", "y.values")) %>%
                       map_dfc(1) %>%
                       ##Get cutoff that maximizes Accuracy
                       filter(y.values==max(y.values)) %>%
                       pull(x.values) %>%
                       median() %>%
                       tibble(BestCutoff = .))

      ##Use cutoffs for predictions
      pred <- pred %>%
        left_join(cutoffs, by=cutoffBy) %>%
        mutate(Predicted = if_else(ProbOfProbClass >= BestCutoff,
                                   probClass,
                                   nonProbClass)) %>%
        select(rowIndex, Actual, Predicted,
               ##Rename ProbOfProbClass
               "Prob{probClass}" := ProbOfProbClass, {{group}},
               !!cutoffBy, BestCutoff)

      ##Add pred class attribute
      class(pred) <- c("pred", class(pred))

      return(pred)
    }
  }


  # Check arguments ---------------------------------------------------------
  ##x can either be a classifier or a predictions dataframe
  if (!any(c("train", "pred") %in% class(x))) {
    stop("x must be either a classifier run with caret::train(), ",
         "or a predictions dataframe created by cls_fairness(x, output='pred')")
  }

  ##Match arg options
  output <- match.arg(output)
  cutoffBy <- match.arg(cutoffBy)


  # Generate output -----------------------------------------------------------
  ##Calculate (if needed) and return (optionally) predictions dataframe
  if ("pred" %in% class(x)) {
    if (output=="pred") {
      warning("output=\"pred\" but x is already a predictions dataframe\nReturning x")
      return(x)
    }
    pred <- x
  } else {
    pred <- getPred(x)
  }

  ##If output is pred, return predictions df
  if (output=="pred") {
    return(pred)
  }

  ##If output is dataframe, get df
  if (output=="dataframe") {
    suppressWarnings(suppressMessages(library(tidyr)))

    ##Optionally calculate accuracy by class
    if (!byClass) {
      df <-
        pred %>%
        mutate(Correct = Actual==Predicted) %>%
        group_by({{group}}) %>%
        summarise(across(Correct, mean)) %>%
        pivot_wider(names_from={{group}}, values_from=Correct)
    } else {
      df <-
        pred %>%
        group_by({{group}}, Actual) %>%
        summarise(ClassAcc = mean(Actual==Predicted), .groups="drop") %>%
        pivot_wider(names_from={{group}}, values_from=ClassAcc)
    }

    return(df)
  }

  ##If output is chisq, calculate chisq
  if (output=="chisq") {
    if (!byClass) {
      ##Get hit/miss matrix
      hitMiss <-
        pred %>%
        mutate(Correct = Actual==Predicted) %>%
        group_by({{group}}) %>%
        summarise(Hit = sum(Correct), Miss = sum(!Correct)) %>%
        select(-{{group}}) %>%
        as.matrix()

      ##Optionally rescale token counts to adjust for resampling
      predictionsPerToken <- nrow(pred) / n_distinct(pred$rowIndex)
      if (unResample) {
        hitMiss <- hitMiss / predictionsPerToken
        if (predictionsPerToken==1) {
          message("Predictions dataframe contains just one prediction per ",
                  "token.\nIgnoring unResample=TRUE.")
        }
      } else if (predictionsPerToken > 1) {
        warning("Predictions dataframe contains ",
                round(predictionsPerToken, 1),
                " predictions per token due to resampling.\n",
                "Using output=\"chisq\" with unResample=FALSE can lead to ",
                "Type I errors.")
      }

      ##Calculate and return chisq test
      chisq <- prop.test(hitMiss)
    } else {
      suppressWarnings(suppressMessages(library(purrr)))

      ##Get list of hit/miss matrices
      hitMiss <-
        pred %>%
        count(Actual, Predicted, {{group}}) %>%
        mutate(Type = if_else(Actual==Predicted, "Hit", "Miss")) %>%
        select(-Predicted) %>%
        pivot_wider(names_from=Type, values_from=n) %>%
        select(-{{group}}) %>%
        arrange(Actual) %>%
        nest(data = c(Hit, Miss)) %>%
        pull(data, Actual) %>%
        map(as.matrix)

      ##Optionally rescale token counts to adjust for resampling
      predictionsPerToken <- nrow(pred) / n_distinct(pred$rowIndex)
      if (unResample) {
        hitMiss <- hitMiss %>% map(~ .x / predictionsPerToken)
        if (predictionsPerToken==1) {
          message("Predictions dataframe contains just one prediction per ",
                  "token.\nIgnoring unResample=TRUE.")
        }
      } else if (predictionsPerToken > 1) {
        warning("Predictions dataframe contains ",
                round(predictionsPerToken, 1),
                " predictions per token due to resampling.\n",
                "Using output=\"chisq\" with unResample=FALSE can lead to ",
                "Type I errors.")
      }

      ##Calculate and return chisq tests
      chisq <- map(hitMiss, prop.test)
    }

    return(chisq)
  }

  ##If output is cm, return confusion matrix
  if (output=="cm") {
    suppressWarnings(suppressMessages(library(purrr)))

    cm <-
      pred %>%
      select({{group}}, Predicted, Actual) %>%
      nest(data = -{{group}}) %>%
      arrange({{group}}) %>%
      pull(data, {{group}}) %>%
      map(table)

    ##Optionally rescale token counts to adjust for resampling
    if (unResample) {
      predictionsPerToken <- nrow(pred) / n_distinct(pred$rowIndex)
      cm <- cm %>%
        map(~ .x / predictionsPerToken)
    }

    ##Return cm
    return(cm)
  }
}

##Generate one-row dataframe of fairness & (optionally) performance metrics
cls_summary <- function(x,
                        ##performance (logical): Include columns for non-grouped
                        ##  performance metrics? If TRUE (default), x cannot be
                        ##  predictions dataframe
                        performance=TRUE,
                        group=race,
                        ##refGroup (character): Value of group variable to serve
                        ##  as "reference group". Reference group columns are to
                        ##  the left in output dataframe, and difference columns
                        ##  are calculated as refGroup minus non-refGroup
                        refGroup="black",
                        cutoffBy=c("Resample","race"),
                        unResample=TRUE) {
  suppressWarnings(suppressMessages(library(dplyr)))
  suppressWarnings(suppressMessages(library(tidyr)))
  suppressWarnings(suppressMessages(library(purrr)))

  ##Check x class
  if (!any(c("train", "pred") %in% class(x))) {
    stop("x must be either a classifier run with caret::train(), ",
         "or a predictions dataframe created by cls_fairness(x, output='pred')")
  }

  ##x can be a predictions dataframe if performance=FALSE
  if ("pred" %in% class(x)) {
    if (performance) {
      stop("If performance=TRUE, x must be a classifier run with caret::train()")
    }
    pred <- x
  }

  ##If x is a classifier, get pred now to save time
  if ("train" %in% class(x)) {
    pred <- cls_fairness(x, "pred", group={{group}}, cutoffBy=cutoffBy)
  }

  ##Check that refGroup is one of the groups
  groups <-
    pred %>%
    pull({{group}}) %>%
    unique()
  if (!(refGroup %in% groups)) {
    stop("refGroup (", refGroup, ") is not one of the groups in ",
         as_name(enquo(group)), ": ", paste(groups, collapse=" & "))
  }

  smry <-
    cbind(
      ##Overall accuracy
      cls_fairness(pred, "dataframe", group={{group}}) %>%
        relocate(all_of(refGroup), .before=1) %>%
        mutate(Diff = .[[1]] - .[[2]]) %>%
        rename_with(~ paste0("Acc_", .x)),

      ##Class accuracies
      cls_fairness(pred, "dataframe", byClass=TRUE, group={{group}}) %>%
        relocate(all_of(refGroup), .after=1) %>%
        mutate(Diff = .[[2]] - .[[3]]) %>%
        pivot_wider(names_from=Actual, names_glue="ClassAcc_{Actual}_{.value}",
                    values_from=-Actual),

      ##Overall accuracy chisq
      cls_fairness(pred, "chisq", unResample=unResample) %>%
        with(tibble(Acc_Chisq_stat = statistic, Acc_Chisq_df = parameter,
                    Acc_Chisq_p = p.value)),

      ##Class accuracy chisq
      cls_fairness(pred, "chisq", byClass=TRUE, unResample=unResample) %>%
        map_dfr(~ with(.x, tibble(Chisq_stat = statistic,
                                  Chisq_df = parameter,
                                  Chisq_p = p.value)),
                .id="Class") %>%
        pivot_wider(names_from=Class, names_glue="ClassAcc_{Class}_{.value}",
                    values_from=-Class))

  ##Optionally add performance
  if (performance) {
    smry <- smry %>%
      cbind(x %>%
              pluck("resample") %>%
              summarise(across(matches("Acc|AUC"), mean)) %>%
              rename_with(~ gsub("Accuracy", "Acc", .x)))
  }

  smry
}


# Display utilities -----------------------------------------------------------

##Print data
printData <- function(data,
                      strategy,
                      dependent=Rpresent,
                      group=race,
                      outlierCols=FALSE,
                      excludeCols=c(Rpresent, race, contains("Outlier")),
                      umsTxt="../Input-Data/UMS-List.txt") {
  ##Validate args
  umsValidate(strategy, umsTxt=umsTxt)
  stopIfNoCol(data, {{dependent}}, "dependent", "data")
  stopIfNoCol(data, {{group}}, "group", "data")
  stopIfNoCol(data, {{excludeCols}}, "excludeCols", "data")

  ##Print counts
  suppressWarnings(suppressMessages(library(dplyr)))
  dataSmall <- data %>%
    select(-{{excludeCols}})
  tableDepGp <- data %>%
    select({{dependent}}, {{group}}) %>%
    table()
  header <- paste("Data:", nrow(dataSmall), "tokens with",
                  ncol(dataSmall), "predictors")
  if (outlierCols) {
    numOutliers <- sum(grepl("Outlier", colnames(data)))
    header <- paste(header, "and", numOutliers, "outlier columns")
  }
  cat(header, fill=TRUE)
  print(addmargins(tableDepGp))

  ##Add additional proportion display for proportional UMSs (or baseline, for
  ##  comparison's sake)
  if (strategy %in% c("0.0", "1.3.1", "1.3.2", "1.4") ||
      startsWith(strategy, "4")) {
    gpname <- as_name(enquo(group))
    cat("Proportionally by ", gpname, ":", sep="", fill=TRUE)
    print(addmargins(prop.table(tableDepGp, 2), 1))
  }
}

##Print generic header (or separator)
printHeader <- function(text, char="-", width=80, pre="\n", post="\n") {
  cat(pre, paste(rep(char, width), collapse=""), "\n", sep="", fill=FALSE)
  if (!missing(text)) {
    cat(text, fill=TRUE)
    cat(paste(rep(char, width), collapse=""), post, sep="", fill=FALSE)
  }
}

##Print UMS header
printHeaderUMS <- function(strategy, umsTxt="../Input-Data/UMS-List.txt", ...) {
  ##Validate strategy & get UMS list
  umsList <- umsValidate(strategy, umsTxt=umsTxt)
  description <- umsList$Description[umsList$UMS==strategy]

  ##Handle combination UMSs
  if (startsWith(description, "Combination")) {
    substrat1 <- sub("Combination of (.+) & (.+)", "\\1", description)
    substrat2 <- sub("Combination of (.+) & (.+)", "\\2", description)
    subdesc1 <- umsList$Description[umsList$UMS==substrat1]
    subdesc2 <- umsList$Description[umsList$UMS==substrat2]
    description <- paste0("Combination of:",
                          "\n  - ", substrat1, ": ", subdesc1,
                          "\n  - ", substrat2, ": ", subdesc2)
  }

  ##Put together header (pass ... args to printHeader())
  hTxt <- paste0("Unfairness mitigation strategy ", strategy, ":\n  ",
                 description)
  printHeader(hTxt, ...)
}


# Meta-utilities --------------------------------------------------------------

##Check that umsTxt path exists in working dir (using script dir as fallback)
umsPath <- function(umsTxt) {
  suppressWarnings(suppressMessages(library(this.path)))
  scriptDir <- this.dir()
  relPath <- path.join(scriptDir, umsTxt)
  if (file.exists(umsTxt)) {
    return(umsTxt)
  } else if (file.exists(relPath)) {
    return(relPath)
  } else {
    stop("Could not find umsTxt in...\n",
         "* working directory (", getwd(), ")\n",
         "* directory containing UMS-Utils.R (", scriptDir, ").\n",
         "Change umsTxt or call from different working directory.")
  }
}

##Check that strategy is valid
umsValidate <- function(strategy, umsTxt) {
  ##Check strategy against valid UMSs
  umsList <- read.csv(umsPath(umsTxt), sep="\t")
  validStrategy <- umsList$UMS
  if (!(strategy %in% validStrategy)) {
    stop("strategy ", strategy, " not recognized")
  }

  ##"precursor models" shouldn't be optimized for performance
  if (grepl("0\\.[^0]", strategy)) {
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    fileArg <- "--file="
    scriptName <- sub(fileArg, "", cmdArgs[grep(fileArg, cmdArgs)])
    if (length(scriptName) != 0 &&
        scriptName %in% c("Hyperparam-Tuning.R", "Outlier-Dropping.R")) {
      stop(scriptName, " can't be used for precursor strategy ", strategy)
    }
  }

  ##Return umsList invisibly
  invisible(umsList)
}

##Helper function to check existence of columns (using tidyselect semantics)
stopIfNoCol <- function(data, checkCol, colLabel, datLabel="training data") {
  suppressWarnings(suppressMessages(library(dplyr)))
  suppressWarnings(suppressMessages(library(rlang)))
  tryCatch(data %>% select({{checkCol}}),
           error = function(e) stop(colLabel, " argument (",
                                    as_name(enquo(checkCol)),
                                    ") does not exist in ", datLabel))
}

##Helper function for print-debugging
announce <- function(x) {
  nm <- deparse(substitute(x))
  cat(nm, fill=TRUE)
  print(x)
  cat(fill=TRUE)
}



##Display helptext if this script called directly
if ("--file=UMS-Utils.R" %in% commandArgs()) {
  cat("UMS-Utils.R:",
      paste0("   ", c(
        "Utility functions for generating and analyzing auto-coders that utilize",
        "unfairness mitigation strategies.",
        "\n",
        "Not meant to be run directly with Rscript, but can be source()d so utility",
        "functions can be used interactively."
      )), fill=TRUE)
}
