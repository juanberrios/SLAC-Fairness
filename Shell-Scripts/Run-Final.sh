#!/usr/bin/env bash

#SBATCH --job-name=Run-Final
#SBATCH --qos=short
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=3
#SBATCH --output=../Outputs/Shell-Scripts/Run-Final.out

##EDITABLE PARAMETERS:
##Baseline UMS
UMS_base=0.0
##Fair UMS
UMS_fair=4.2.1
##Optimal tuning grid: baseline UMS
tgrid_base="data.frame(mtry = 15, splitrule = 'gini', min.node.size = 1)"
##Optimal tuning grid: fair UMS
tgrid_fair="data.frame(mtry = 15, splitrule = 'gini', min.node.size = 5)"
##Outliers to drop: baseline UMS
drop_base=F2max_Outlier,F1_55_Outlier,F3_45_Outlier,F1_80_Outlier,F3_40_Outlier
##Outliers to drop: fair UMS
drop_fair=F1_55_Outlier,F3_80_Outlier,F3_75_Outlier,F2_45_Outlier,F2_80_Outlier,F1_25_Outlier
##Path to directory containing R script(s)
scriptdir=../R-Scripts/
##Path *from scriptdir* to tmpfile storing package list
ppath=../Outputs/Other/packages_final.tmp 
##Path *from scriptdir* to directory for saving final classifiers
savedir=../Outputs/Autocoders-to-Keep/

##Header
sep="\n==============================================================================\n"
printf "$sep"
printf "Final classifiers:\n"
printf "Baseline:\n"
printf "  UMS: $UMS_base\n"
printf "  Tuning grid: $tgrid_base\n"
printf "  Outliers dropped: $drop_base\n"
printf "Fair:\n"
printf "  UMS: $UMS_fair\n"
printf "  Tuning grid: $tgrid_fair\n"
printf "  Outliers dropped: $drop_fair"
printf "$sep"

##Load modules (if using Lmod)
if command -v module &> /dev/null ; then module load gcc/10.2.0 r/4.2.0 ; fi

##Run files
cd $scriptdir
Rscript Run-UMS.R --ums $UMS_base --drop "$drop_base" --tune-grid "$tgrid_base" --output-prefix Final-Cls --output-dir $savedir --pkglist-path $ppath
Rscript Run-UMS.R --ums $UMS_fair --drop "$drop_fair" --tune-grid "$tgrid_fair" --output-prefix Final-Cls --output-dir $savedir --pkglist-path $ppath

##Finish up
Rscript Session-Info.R --pkglist-path $ppath
rm $ppath
##Run job stats script if it exists
command -v crc-job-stats &> /dev/null && command crc-job-stats
