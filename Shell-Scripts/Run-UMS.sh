#!/usr/bin/env bash

#SBATCH --job-name=Run-UMS
#SBATCH --qos=short
#SBATCH --time=04:00:00
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=3
#SBATCH --output=../Outputs/Shell-Scripts/Run-UMS.out

##Parse args (first is UMS code, second is flag to optimize for performance)
UMS=${1}
optimize=${2}

##EDITABLE PARAMETERS:
##Path to directory containing R script(s)
scriptdir=../R-Scripts/
##Path *from scriptdir* to tempfile storing package list
ppath=../Outputs/Other/packages_UMS$UMS.tmp
#### Only apply if optimizing for performance (flag -o) ####
##Path *from scriptdir* to csv storing optimal tuning parameters
tgrid=../Outputs/Other/Best-Params_UMS$UMS.csv
##Path *from scriptdir* to tempfile storing dropped outliers
dropped=../Outputs/Other/dropped-outliers_UMS$UMS.tmp
##Path *from scriptdir* to directory for saving optimized classifier (set to NULL to skip saving)
savedir=../Outputs/Diagnostic-Files/Temp-Autocoders/
##File prefix for saving optimized classifier
saveprefix=Run-Cls

##Check args
if [ -z $UMS ]; then
  echo You must specify a UMS
  echo "Syntax: bash Run-UMS.sh UMS [-o]"
  exit 1
fi
if [ $UMS = "-o" ]; then
  echo You must specify a UMS
  echo "Syntax: bash Run-UMS.sh UMS [-o]"
  exit 1
fi
if [ -n "${optimize}" ] && [ $optimize != "-o" ]; then
  echo Second option can only be -o
  echo "Syntax: bash Run-UMS.sh UMS [-o]"
  exit 1
fi

##Load modules (if using Lmod)
if command -v module &> /dev/null ; then module load gcc/10.2.0 r/4.2.0 ; fi

##Run auto-coder
cd $scriptdir
Rscript Run-UMS.R --ums $UMS --pkglist-path $ppath

##Optionally optimize for performance
if [ -n "${optimize}" ]; then
  ##Optimize
  Rscript Hyperparam-Tuning.R --ums $UMS --best-tune-grid $tgrid --pkglist-path $ppath
  Rscript Outlier-Dropping.R --ums $UMS --tune-grid $tgrid --best-drop-cols $dropped --pkglist-path $ppath
  ##Run final auto-coder
  Rscript Run-UMS.R --ums $UMS --tune-grid $tgrid --drop-cols $dropped --output-dir $savedir --output-prefix $saveprefix --pkglist-path $ppath
fi

##Finish up
Rscript Session-Info.R --pkglist-path $ppath

rm $ppath
##Run job stats script if it exists
command -v crc-job-stats &> /dev/null && command crc-job-stats
