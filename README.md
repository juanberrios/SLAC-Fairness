# `SLAC-Fairness`: Tools to assess fairness and mitigate unfairness in sociolinguistic auto-coding

_Dan Villarreal (Department of Linguistics, University of Pittsburgh)_

![](https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png)

This work is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-nc-sa/4.0/).


## Introduction

This GitHub repository is a companion to the paper "Sociolinguistic auto-coding has fairness problems too: Measuring and mitigating overlearning bias", forthcoming in _Linguistics Vanguard_.
In the paper, I investigate **sociolinguistic auto-coding (SLAC)** through the lens of **machine-learning fairness**.
Just as some algorithms produce biased predictions by _overlearning_ group characteristics, I find that the same is true for SLAC.
As a result, I attempt **unfairness mitigation strategies (UMSs)** as techniques for removing gender bias in auto-coding predictions (without harming overall auto-coding performance too badly).


**_Repository navigation:_**

- [_Repository homepage_](https://djvill.github.io/SLAC-Fairness)
- [_Repository code_](https://github.com/djvill/SLAC-Fairness)
- [_Analysis walkthrough_](https://djvill.github.io/SLAC-Fairness/Analysis-Walkthrough)
- [_Unfairness mitigation strategy descriptions_](https://djvill.github.io/SLAC-Fairness/UMS-Info)


### If you're new to sociolinguistic auto-coding (SLAC)

Sociolinguistic auto-coding is a machine-learning method for classifying variable linguistic data (often phonological data), such as the alternation between _park_ & "_pahk_" or _working_ & _workin'_.

You can learn more about SLAC by reading the following resources.

- Dan Villarreal et al.'s 2020 _Laboratory Phonology_ article ["From categories to gradience: Auto-coding sociophonetic variation with random forests"](https://doi.org/10.5334/labphon.216)
- Tyler Kendall et al.'s 2021 _Frontiers in AI_ article ["Considering performance in the automated and manual coding of sociolinguistic variables: Lessons from variable (ING)"](https://doi.org/10.3389/frai.2021.648543)
- Dan Villarreal et al.'s 2019 tutorial ["How to train your classifier"](https://nzilbb.github.io/How-to-Train-Your-Classifier/How_to_Train_Your_Classifier.html)
  - Explains the R code powering my implementation of SLAC


## What's the point of this repository?

First, you can **reproduce** the analysis I performed for the _Linguistics Vanguard_ paper, using the same data and code that I did.
Simply follow the analysis walkthrough [tutorial](https://djvill.github.io/SLAC-Fairness/Analysis-Walkthrough.html).

Second, you can also [adapt this code](#adapting-this-code-to-your-own-projects) to your own projects.
You might want to use it if you want to (1) **assess fairness** for a [pre-existing auto-coder](#assessing-fairness-for-a-pre-existing-auto-coder) and/or (2) create a **fair auto-coder** by [testing unfairness mitigation strategies](#testing-unfairness-mitigation-strategies) on your data.

Finally, I invite [comments, critiques, and questions](#auditing-this-code-to-critique-andor-suggest-changes) about this code.
I've made this code available for transparency's sake, so please don't hesitate to reach out!


## What's in this repository?

The files in this repository fall into a few categories.
Click the links below to jump to the relevant subsection:

- Info written for humans
  - `README.md`: What you're reading now
  - `Analysis-Walkthrough.Rmd` & `.html`: Tutorial for replicating analysis in paper
  - `UMS-Info.Rmd` & `.html`: Descriptions of unfairness mitigation strategies
- [Input data](#input-data)
  - `Input-Data/`
- [Outputs](#outputs)
  - `Outputs/`
- [Code that does stuff](#code)
  - `R-Scripts/`
  - `Shell-Scripts/`
- [Code/info pertaining to the repository itself](#codeinfo-pertaining-to-the-repository-itself)
  - `.gitignore`
  - `LICENSE.md`
  - `renv/`
  - `renv.lock`
  - `.Rprofile`

You can browse files [here](https://github.com/djvill/SLAC-Fairness).


### A quick note on the two-computer setup

This repository's structure reflects the two-computer setup I used to run this analysis.
I generated and measured auto-coders on a more powerful system that is not quite as user-friendly (Pitt's [CRC](https://crc.pitt.edu/)), then analyzed the metrics on my less-powerful-but-user-friendlier laptop.
(It's perfectly fine to use a one-computer setup if you don't have access to high-performance computing;
the code [will take longer to run](https://djvill.github.io/SLAC-Fairness/Analysis-Walkthrough.html#running-time) on a less-powerful machine, but it still might be faster/easier than the HPC learning curve!)
In the rest of this section, I'll refer to this two-computer split several times.


### Input data

Contents:

- `Input-Data/`
  - `LabPhonClassifier.Rds`: Pre-existing auto-coder to analyze for fairness. This auto-coder is the same as in ["How to train"](https://github.com/nzilbb/How-to-Train-Your-Classifier/blob/main/LabPhonClassifier.Rds), but with a `Gender` column added to the auto-coder's `trainingData` element.
  - `trainingData.Rds`: /r/ data for generating auto-coders that use unfairness mitigation strategies, also available [here](https://github.com/nzilbb/How-to-Train-Your-Classifier/blob/crc-version/Data/trainingData.Rds). This is the result of [step 1 in "How to train"](https://nzilbb.github.io/How-to-Train-Your-Classifier/How_to_Train_Your_Classifier.html#step-1). A dataset with more tokens (but less acoustic information) is available [here](https://github.com/nzilbb/Sld-R-Data).
  - `meanPitches.csv`: Pitch data to use for UMS 3.1 (normalizing speaker pitch). I measured pitch (F0) for word-initial /r/ tokens and calculated each speaker's average minimum and maximum pitch.
  - `UMS-List.txt`: Tab-separated file matching UMS codes to descriptions. This is also used by the R scripts to define the set of acceptable UMS codes.

The /r/ and pitch data comes from Southland New Zealand English, historically New Zealand's only regional variety, which is characterized by variable rhoticity.
The New Zealand Institute of Language, Brain and Behaviour maintains a corpus of sociolinguistic interviews with Southland English speakers totaling over 83 hours of data.
This corpus is hosted in an instance of [LaBB-CAT](https://labbcat.canterbury.ac.nz/);
the data files were downloaded from LaBB-CAT, with subsequent data-wrangling in R (including speaker anonymization).

_Skip ahead for info on using your own [auto-coder](#assessing-fairness-for-a-pre-existing-auto-coder) and [training data](#using-your-own-training-data), or modifying the set of [UMSs](#adding-andor-subtracting-umss)._


### Outputs

Contents:

- `Outputs/`
  - `Autocoders-to-Keep/`
    - "Final" auto-coders (saved as `.Rds` files). Unlike the temporary auto-coders, this folder is version-controlled (see [info on `.gitignore`](#codeinfo-pertaining-to-the-repository-itself)), so it's useful for selectively saving auto-coders we want to hold onto.
  - `Shell-Scripts/`
    - Text files (saved with the `.out` file extension) that record any output of [shell scripts](#shell-scripts), including errors. Useful for diagnosing issues with the code if something goes wrong.
  - `Performance/`
    - Tabular data (saved as `.csv` files) with metrics of auto-coders' performance (e.g., overall accuracy) and fairness (e.g., accuracy for women's vs. men's tokens). These files bridge the [two-computer split](#a-quick-note-on-the-two-computer-setup) split: we extract metrics on a more powerful system (see [walkthrough](https://djvill.github.io/SLAC-Fairness/Analysis-Walkthrough.html#baseline-metrics)) so we can analyze them on a user-friendlier computer.
  - `Diagnostic-Files/`: Temporary files that are useful only in the moment (e.g., peeking "under the hood" to diagnose a code issue if something goes wrong) and/or too large to [share between computers](#a-quick-note-on-the-two-computer-setup). Most files are [`.gitignore`d](#codeinfo-pertaining-to-the-repository-itself), save for empty `dummy_file`s that exist only so the [empty folders can be shared to GitHub](https://stackoverflow.com/a/8418403).
    - `Model-Status/`
      - Temporary files (with extension `.tmp`) that are created during [optimization for performance](https://djvill.github.io/SLAC-Fairness/Analysis-Walkthrough.html#baseline) to signal which auto-coders are completed or running.
    - `Temp-Autocoders/`
      - Auto-coders run with different UMSs, for which we want to measure performance and fairness but we don't need to version-control
  - `Other/`: Files mostly meant for passing info between scripts
    - `Var-Imp*.csv`: Data on variable importance for ["precursor" UMSs](https://djvill.github.io/SLAC-Fairness/Analysis-Walkthrough.html#meas-precursor) used for UMSs 2.1.x.
    - `Best-Params*.csv`: Optimal hyperparameters from [hyperparameter-tuning](https://nzilbb.github.io/How-to-Train-Your-Classifier/How_to_Train_Your_Classifier.html#step-4) runs.
    - `Drop-Log*.csv`: Log files for [outlier-dropping](https://nzilbb.github.io/How-to-Train-Your-Classifier/How_to_Train_Your_Classifier.html#step-5) runs.

All these outputs were generated using the data and code in this repository.
You may want to create a 'clean' version of the repository without any of these outputs, to see if your system replicates the outputs I got.


### Code

Contents:

- `R-Scripts/`: Scripts that do the heavy lifting of running auto-coders and facilitating analysis.
- `Shell-Scripts/`: Scripts meant for the user to run; these scripts call the R scripts and collate their outputs.


The division of labor has two benefits.
First, it makes the code more modular, so a larger process isn't completely lost if just one part fails.
Second, many high-performance computing environments don't allow users to run code on-demand, instead submitting job requests, packaged into shell scripts, to a workload management system (aka job queue).
These shell scripts are written to be compatible with [Slurm](https://slurm.schedmd.com/), the job queue used by Pitt's [CRC](https://crc.pitt.edu/) clusters.
If your computing environment _doesn't_ require submitting job requests, the shell scripts should still run as-is.
You also have the option of foregoing the shell scripts and running the R scripts directly.




#### R scripts

These include 'main scripts' that run auto-coders and 'helper scripts' that define functionality shared among the main scripts.

Main scripts:

- `R-Scripts/`
  - `Run-UMS.R`: Generates a single auto-coder according to an unfairness mitigation strategy.
  - `Hyperparam-Tuning.R`: Subjects an auto-coder to [hyperparameter tuning](https://nzilbb.github.io/How-to-Train-Your-Classifier/How_to_Train_Your_Classifier.html#step-4), one stage of optimizing an auto-coder for performance. (Note: This tunes only what "How to train" calls [`ranger` parameters](https://nzilbb.github.io/How-to-Train-Your-Classifier/How_to_Train_Your_Classifier.html#ranger-parameters-mtry-splitrule-min.node.size)) because I'm now less sure that the other hyperparameters are appropriate for tuning.)
  - `Outlier-Dropping.R`: Subjects an auto-coder to [outlier dropping](https://nzilbb.github.io/How-to-Train-Your-Classifier/How_to_Train_Your_Classifier.html#step-5), one stage of optimizing an auto-coder for performance.

The main scripts are written to be called from a command-line client like Bash, using the command `Rscript`.
(To use `Rscript`, R needs to be in your [PATH](#running-this-code-on-your-own-machine).)
For example, if you navigate Bash to `R-Scripts/`, you can run `Rscript Run-UMS.R --ums 0.0`
These scripts take several arguments (like `--ums`);
to see arguments, run `Rscript <script-name> --help` from the command line.
If you prefer working exclusively in R, you can use `rscript()` from [`callr`](https://cran.r-project.org/package=callr) to call these scripts from within R (e.g., `callr::rscript("Run-UMS.R", c("--ums", "0.0"), wd="R-Scripts/", stdout="../Shell-Scripts/Output/Run-UMS_UMS0.0.out", stderr="2>&1")`).


Helper scripts:

- `R-Scripts/`
  - `UMS-Utils.R`: Contains utility functions for generating and analyzing auto-coders that utilize UMSs. The most important functions are:
    - `umsData()`: Reshape data for auto-coder by applying UMS, only keeping necessary columns, and optionally dropping outliers
    - `umsFormula()`: Specify model formula based on UMS
    - `cls_fairness()`: Investigate auto-coder fairness (see [walkthrough](#rq2-ums-utils))
    - `cls_summary()`: Generate one-row dataframe of fairness/performance metrics (see [walkthrough](#rq2-ums-utils))
  - `Rscript-Opts.R`: Defines command-line options for how main scripts should run.
  - `Session-Info.R`: Combines and prints R session info from the outputs of multiple scripts. Meant to be used in shell scripts.


#### Shell scripts

Contents:

- `Shell-Scripts/`
  - `Run-UMS.sh`: Generates a single auto-coder according to a UMS, and optionally optimizes it for performance. This flexible lower-level script is useful for exploratory analysis.
  - `Baseline.sh`: Wrapper script that calls `Run-UMS.sh` for [baseline auto-coder](https://djvill.github.io/SLAC-Fairness/Analysis-Walkthrough.html#baseline) (mostly exists to override the default [outfile](#outputs) name)
  - `UMS-Round1.sh`: Generates auto-coders according to UMSs whose codes start with 0, 1, 2, or 3 (save for UMS 0.0, the baseline).
  - `UMS-Round2.sh`: Generates auto-coders according to UMSs whose codes start with 4.

The shell scripts are written to be called from Bash, using the commands `bash` (to run directly) or `sbatch` (to submit to a Slurm job queue; [see above](#code)).
For example, if you navigate Bash to `Shell-Scripts/`, you can run `sbatch Baseline.sh`.
`Run-UMS.sh` takes two arguments: a UMS numerical code (e.g., `sbatch Run-UMS.sh 4.2.1`) and an optional `-o` flag to [optimize the auto-coder for performance](https://djvill.github.io/SLAC-Fairness/Analysis-Walkthrough.html#baseline).
All other shell scripts hard-code these options (as well as other options passed to the R scripts), but these can be adjusted as needed (under the heading `##EDITABLE PARAMETERS`).


The shell scripts also collate output from the R scripts;
I recommend saving this output to a text file.
These scripts include a Slurm command that automatically writes script output to a corresponding `.out` file in `Outputs/Shell-Scripts/`.
If you're not using Slurm, you can append a command that tells Bash where to send outputs, including errors (e.g., `bash Baseline.sh &> ../Outputs/Shell-Scripts/Baseline.out`);
if you omit this part of the command (e.g., `bash Baseline.sh`), the output will simply print in Bash.


CRC's cluster uses [Lmod](http://lmod.readthedocs.org) to make modules (like R) available to shell scripts via the `module load` command.
Your system may not need to load modules explicitly, or may use different commands to load R.


### Code/info pertaining to the repository itself

Contents:

- `.gitignore`: Tells Git which files/folders to exclude from being version-controlled (and being shared to GitHub or [between computers](#a-quick-note-on-the-two-computer-setup)). Because the auto-coders are huge files, I exclude `Outputs/Diagnostic-Files/Temp-Autocoders/` from version-control and just [pull out fairness/performance data](https://djvill.github.io/SLAC-Fairness/Analysis-Walkthrough.html#baseline-measure) instead. If there's any I want to keep, I save them to the non-ignored folder `Outputs/Autocoders-to-Keep/`.
- `LICENSE.md`: Tells you what you're permitted to do with this code.
- `renv/`: Set up by the [`renv` package](https://rstudio.github.io/renv/) to ensure our code behaves the same regardless of package updates. See more info [below](#renv).
- `renv.lock`: Set up by `renv` to store [info about package versions](https://rstudio.github.io/renv/articles/lockfile.html).
- `.Rprofile`: Contains R code to run at the start of any R session in this repository. In this case, this code was set up by `renv` to run a script that loads the package versions recorded in `renv.lock`. If you want to disable `renv`, simply delete this file.



## Running this code on your own machine

To run this code on your own machine, you'll need a suitable computing environment and software.
All required and recommended software is free and open-source.
This document was originally run using high-performance computing resources provided by the University of Pittsburgh's [Center for Research Computing (CRC)](https://crc.pitt.edu/), in particular its [shared memory parallel cluster](https://crc.pitt.edu/resources/h2p-user-guide/node-configuration).
You _can_ run this code on a normal desktop or laptop---it just might take a while!
You'll also need at least 400 Mb of disk space free.
See the [walkthrough](https://djvill.github.io/SLAC-Fairness/Analysis-Walkthrough.html#script-info) for more information about machine specs, running time, and disk space used.


Required software:

- The statistical computing language [R](https://cloud.r-project.org/) (version >= 4.3.1)
  - Since these scripts call R from the command line, R must be in your PATH (directions for [Windows](https://docs.oracle.com/en/database/oracle/machine-learning/oml4r/1.5.1/oread/creating-and-modifying-environment-variables-on-windows.html#GUID-DD6F9982-60D5-48F6-8270-A27EC53807D0), [macOS](https://www.architectryan.com/2012/10/02/add-to-the-path-on-mac-os-x-mountain-lion/#.Uydjga1dXDg), [Unix](https://unix.stackexchange.com/a/26059)). If R is in your PATH, then at the command line `Rscript -e R.version.string` will print your R version.
- R packages:
  - `tidyverse` (v. >= 2.0.0)
  - `magrittr` (v. >= 2.0.3)
  - `caret` (v. >= 6.0-94)
  - `ranger` (v. >= 0.15.1)
  - `ROCR` (v. >= 1.0-11)
  - `foreach` (v. >= 1.5.2)
  - `doParallel` (v. >= 1.0.17)
  - `optparse` (v. >= 1.7.3)
  - `this.path` (v. >= 1.4.0.13, install via `remotes::install_github("ArcadeAntics/this.path", "5d755e1")`)
  - `benchmarkme` (v. >= 1.0.8)
  - `lubridate` (v. >= 1.9.2)
  - `rmarkdown` (v. >= 2.22)
  - `knitr` (v. >= 1.43)
  - `renv` (v. >= 0.17.3)
  - These packages will install dependencies that you don't need to install directly. See full R session info [here](https://djvill.github.io/SLAC-Fairness/Analysis-Walkthrough.html#R-session-info)
- The command-line client [Bash](https://www.gnu.org/software/bash/) (v. >= 5.0.0)
  - If you install Git (recommended), Bash is included in the install
- The document converter [Pandoc](https://pandoc.org/) (v. >= 2.19)
  - If you install RStudio (recommended), Pandoc is included in the install

Please note that R and its packages are continually updated, so in the future the code may not work as expected (or at all!).
If you hit a brick wall, don't hesitate to [reach out](#auditing-this-code-to-critique-andor-suggest-changes)!


I also recommend using Git and GitHub to create your own shareable version of the code; 
doing so will help me effectively troubleshoot any issues you have. 
In particular:

1. [Download Git](https://git-scm.com/downloads) onto your computer.
    - See this [Git tutorial](https://github.com/djvill/LSA2019-Reproducible-Research) if you've never used it before.
1. Sign up for a free [GitHub account](https://github.com/join)
1. [Fork this repository](https://github.com/djvill/SLAC-Fairness/fork) (keep the same repository name), and clone it onto your computer.
1. Test out the code on your own system: Edit the code, create commits, push your commits to your remote fork.
    - You may want to create a 'clean' version of the repository without any of the [generated outputs](#outputs), to see if your system replicates my outputs.
1. [Reach out](#auditing-this-code-to-critique-andor-suggest-changes)!


Finally, I recommend using the integrated development environment [RStudio](https://www.rstudio.com/products/rstudio/download/).
While it doesn't change how the code in this repository works, RStudio makes R code easier to understand, write, edit, and debug.


### `renv`

This repository uses the [`renv`](https://rstudio.github.io/renv/) package to ensure that updates to R packages don't break the code.
In effect, `renv` freezes your environment in time by preserving the package versions the code was originally run on.
This is great from a reproducibility perspective, but it entails some extra machinery before you can run the code.
For all the examples below, you need to load this repo in R or RStudio by setting your working directory somewhere inside the repo.


Before you can run any of this code, run `renv::restore()`.
This will download the packages at the correct versions to an `renv` cache on your system.
Then you should be able to run this code on your machine.


Of course, using old versions of these packages means you won't be able to benefit to any package updates since this repo was published.
If you want to use new package versions, you have to register them with `renv`.
If you're using R 4.3.x (the version used for this code), run `renv::update()`;
if R >= 4.4, run `renv::init()` and select option 2.
To update renv itself, run `renv::upgrade()`.
Of course, the code may not work as expected thanks to changes to the packages it relies on.


If you want to use a package that's not registered with `renv`, use `renv::record()`.


Finally, if you're finding this all too much of a hassle, you can skip using `renv` altogether;
just delete `.Rprofile` and restart R/RStudio.



## Adapting this code to your own projects

How much you want to adapt this code is really up to you.
You might want to 'carbon-copy' this analysis on your own project, but in all likelihood your project will dictate that you make some changes to better fit your project's needs.
Just as ["training an auto-coder is not a one-size-fits-all process"](https://nzilbb.github.io/How-to-Train-Your-Classifier/How_to_Train_Your_Classifier.html#introduction), so too is auto-coding fairness.
For example, if you are confident that your predictor data (e.g., acoustic measures) does not suffer from measurement error, you can skip the time-consuming step of [accounting for outliers](https://nzilbb.github.io/How-to-Train-Your-Classifier/How_to_Train_Your_Classifier.html#do_we_even_need_to_mark_outliers).
In some cases, this code might not necessarily work for your project; 
for example, this code only handles fairness across two groups, and it only handles binary classification (two categories).

Below, you can read about:

- Assessing fairness for a [pre-existing auto-coder](#assessing-fairness-for-a-pre-existing-auto-coder)
- Creating a fair auto-coder by [testing unfairness mitigation strategies](#testing-unfairness-mitigation-strategies)
  - Using your own [training data](#using-your-own-training-data)
  - Adding and/or subtracting [UMSs](#adding-andor-subtracting-umss)


In addition, I strongly recommend making the data you use for this task publicly available if possible, since open data helps advance science (see Villarreal & Collister "Open methods in linguistics", in press for Oxford collection _Decolonizing linguistics_).
However, if you do so, make sure what you share conforms to the ethics/IRB agreement(s) in place when the data was collected (if applicable).

Finally, if there's anything in this code that you can't figure out or isn't working for you, please **don't hesitate to [reach out](#auditing-this-code-to-critique-andor-suggest-changes)!**
Please note that there is no warranty for this code.


### Assessing fairness for a pre-existing auto-coder

This is one possible goal of your analysis, mirroring the _Linguistics Vanguard_ paper's [RQ2](https://djvill.github.io/SLAC-Fairness/Analysis-Walkthrough.html#RQ2).
To analyze your auto-coder, it needs to have been generated by `caret::train()`.
The auto-coder's `trainingData` element also needs a column with group data (e.g., which tokens belong to female vs. male speakers).
If you use the scripts in this repository, that's taken care of for you; 
`umsData()` retains the group column in the training dataframe passed to `train()`, and `umsFormula()` excludes the group column from the predictor set.
However, if you didn't use these scripts to run your auto-coder, you'll need to either manually add the group column to the `trainingData` element, or just re-run your auto-coder using these scripts.


If your auto-coder conforms to these requirements, you can use functions from `R-Scripts/UMS-Utils.R` to analyze fairness.
See the [walkthrough](https://djvill.github.io/SLAC-Fairness/Analysis-Walkthrough.html#rq2-ums-utils) for examples of how to use this code.


### Testing unfairness mitigation strategies

This is the other possible goal of your analysis, mirroring the _Linguistics Vanguard_ paper's [RQ3](https://djvill.github.io/SLAC-Fairness/Analysis-Walkthrough.html#RQ3).


#### Using your own training data

You'll need your own training data (in place of `trainingData.Rds`), and you may need normalization data depending on which UMSs you want to try.

Formatting requirements for training data:

- Tabular data (data stored in rows and columns), saved in a `.csv` or `.Rds` file
- Each row represents a single token of some categorical linguistic variable
- At least some of the tokens have been coded into classes (in `trainingData.Rds`, these are tokens for which the column `HowCoded=="Hand"`)
- Columns needed for auto-coder:
  - 1 column with variant labels for already-coded tokens and blanks/`NA`s for uncoded tokens (`Rpresent` in `trainingData.Rds`)
    - Currently, this code only handles binary classification (two categories, not counting `NA`s)
  - 1 column with the group that you're assessing fairness for (`Gender` in `trainingData.Rds`)
    - Currently, this code only handles two-group fairness
  - Multiple columns that contain predictors that the auto-coder will use for coding (in `trainingData.Rds`, 180 columns from `tokenDur` to `absSlopeF0`, inclusive)
    - Can be any data type


If you want to perform any speaker normalization (either as a preprocessing step or as UMS 3.1), you'll also need:

- In your training data, a `Speaker` column
- An additional data file with normalization baselines (like `meanPitches.csv`) :
  - One row per speaker, with every speaker in your training data
  - A `Speaker` column
  - A column for each baseline measure you want to use for normalization (`MinPitch` in `meanPitches.csv` is used as baseline for the `F0min` measure in the training data, `MaxPitch` for the `F0max` measure)


In addition to those requirements, here are some data formatting recommendations (you don't _have_ to format your data this way, but if not you'll need to tinker with the code some more).
Some of these pertain to the [data preprocessing step in "How to train"](https://nzilbb.github.io/How-to-Train-Your-Classifier/How_to_Train_Your_Classifier.html#step-1):

- If you suspect that your predictors have [measurement error](https://nzilbb.github.io/How-to-Train-Your-Classifier/How_to_Train_Your_Classifier.html#do_we_even_need_to_mark_outliers) and you want to take advantage of the outlier-dropping script, then you need to mark measurement outliers. Outliers for predictors `X` & `Y` (for example) should be marked as `TRUE` in columns `X_Outlier` & `Y_Outlier`.
- The code drops rows that have `NA`s for the dependent variable, group, and predictor columns (and outlier columns, if applicable). Depending on how many missing measurements you have, you might want to consider imputing measurements and/or thinning your predictor set.
- You might want to add a `HowCoded` column to easily separate hand-coded and auto-coded tokens.
- Thanks to [pre-processing](https://nzilbb.github.io/How-to-Train-Your-Classifier/How_to_Train_Your_Classifier.html#step-1), `trainingData.Rds` reflects normalized measurements for formant timepoints but not pitch, so UMS 3.1 involves pitch normalization. You might decide to fold _all_ normalization into pre-processing and skip normalization as a UMS.
- You may want to anonymize your speakers, especially if you choose to make your data open, as in `trainingData.Rds`, but this is strictly optional.


Finally, as a general note, you may have to tweak the [R scripts](#r-scripts) a little bit to accommodate your data.
For example:

- These scripts assume the training data file is an `.Rds` file. If it's a `.csv`, you'll have to tweak the code
- If your columns have different names than the ones in `trainingData.Rds` (e.g., if your dependent variable isn't `Rpresent`), you'll need to find-and-replace column names in the scripts. Alternatively, you can pass column names as arguments to `UMS-Utils.R` functions (e.g., `umsData(myData, dependent=ING, group=Ethnicity)`).
- If you don't have a `HowCoded` column, you'll need to modify the lines of code that refer to that column.
- Depending on the size of your predictor set, you'll want to change the default value of `mtry` (the number of predictors attempted at each split) in `Rscript-Opts.R` and `Hyperparam-Tuning.R`; a typical value is the square root of the number of predictors, rounded down to the nearest integer.
- If you rename or reorganize folders or files, you'll need to change the code to account for that.


#### Adding and/or subtracting UMSs

Depending on your groups, your predictor set, and/or your dependent variable, you might want to add or subtract UMSs.
For example, if you already have equal token counts for women vs. men, UMS 1.1 (downsample men to equalize token counts by gender) wouldn't apply.


To add a new UMS:

1. Pick a new [UMS code](https://djvill.github.io/SLAC-Fairness/UMS-Info.html#ums-codes)
    - Don't use a code that's already been defined (it just creates unnecessary complications)
    - If it's a combination UMS, the code should start with `4`
1. Add the code and description to `Input-Data/UMS-List.txt`
1. Modify `umsData()` in `R-Scripts/UMS-Utils.R`
    - Single UMSs: Add a new `} else if (UMS=="<new-UMS>") {` block to the `implementUMS()` subroutine
    - Combination UMSs: Add code to interpret the second & third digits near the bottom of `umsData()`
    - Note that `umsData()` uses [`tidyselect` semantics](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html) for several arguments (`dependent`, `group`, `predictors`, & `dropCols`). If you're using any of these column names in a `dplyr` function, wrap them in double-braces (e.g., `data %>% select({{dependent}}, {{group}})`); if you need a column name as a string, use the deparse-substitute trick (e.g., `depName <- deparse(substitute(dependent))`)
1. If using a shell script to run multiple UMSs in a single round, edit the script so the UMS code is matched by the `pattern` regex and not by `excl` (e.g., to include UMS 5.1, use `pattern=^[0-35]`)

You only need to subtract a UMS explicitly if you're using a shell script to run multiple UMSs in a single round.
To subtract a UMS, use the `excl` regex to exclude it (e.g., to exclude UMSs 1.4 and 2.2, use `excl="^1.4|2.2"`).
No need to modify `umsData()`, since the code will just skip over that UMS in the chain of `else if {}` statements.

Note that the existing UMS list is actually more general than its descriptions suggest.
For example, UMSs 1.3.1 and 1.3.2 both achieve equal /r/ base rates by gender, by downsampling either women's Absent (1.3.1) or men's Present (1.3.2).
However, `umsData()` actually translates this into "downsample one of the classes from the smaller group" vs. "the bigger group", automatically detecting which class to downsample from which group.
Try plugging your data into `umsData()` to see whether the existing code affects your data the way you expect.


## Auditing this code to critique and/or suggest changes

Readers are more than welcome to critique this code!
While I think much of this code is pretty solid, there are no doubt some bugs here and there, some inefficient code implementations, and/or some tortured data-scientific reasoning.
You can [raise GitHub issues](https://github.com/djvill/SLAC-Fairness/issues), [start discussions](https://github.com/djvill/SLAC-Fairness/discussions), or [send me an email](mailto:d.vill@pitt.edu?subject=[SLAC-Fairness]%20Auditing%20code).

**Please don't be afraid to suggest changes, report bugs, or ask questions---I want this code to be useful for you, and there are no bad questions!**


## Acknowledgements

I would like to thank Chris Bartlett, the Southland Oral History Project (Invercargill City Libraries and Archives), and the speakers for sharing their data and their voices.
Thanks are also due to Lynn Clark, Jen Hay, Kevin Watson, and the New Zealand Institute of Language, Brain and Behaviour for supporting this research.
Valuable feedback was provided by audiences at NWAV 49, the Penn Linguistics Conference, Pitt Computer Science, and the Michigan State SocioLab.
Other resources were provided by a Royal Society of New Zealand Marsden Research Grant (16-UOC-058) and the University of Pittsburgh Center for Research Computing (specifically, the H2P cluster supported by NSF award number OAC-2117681).
Any errors are mine entirely.

