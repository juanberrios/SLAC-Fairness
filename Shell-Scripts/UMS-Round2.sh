#!/usr/bin/env bash

#SBATCH --job-name=UMS-Round2
#SBATCH --qos=short
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=3
#SBATCH --output=../Outputs/Shell-Scripts/UMS-Round2.out

##EDITABLE PARAMETERS:
##Regex (Perl-compatible) defining UMSs to run
pattern=^4
##Regex (Perl-compatible) defining UMSs to exclude
excl=^$
##Path to directory containing R script(s)
scriptdir=../R-Scripts/
##Path *from scriptdir* to tempfile storing package list
ppath=../Outputs/Other/packages_UMS-Round2.tmp

##Get UMS array
cd $scriptdir
umstxt=$(< ../Input-Data/UMS-List.txt)
umslist=$(echo "$umstxt" | tail -n +2 | grep -oP ".+(?=\t)" | grep -P $pattern | grep -vP $excl)
howmany() { echo $#; }
numums=$(howmany $umslist)
if [ $numums -eq 0 ]; then 
  echo No UMSs with pattern $pattern
  exit 1
fi

##Header
sep="\n==============================================================================\n"
printf "$sep"
printf "Unfairness mitigation strategies: round 2\n"
printf "Running $numums strategies"
printf "$sep"

##Load modules (if using Lmod)
if command -v module &> /dev/null ; then module load gcc/12.2.0 r/4.3.0 ; fi

##Run files (loop over umslist)
for ums in $umslist ; do 
  Rscript Run-UMS.R --ums $ums --pkglist-path $ppath
done

##Finish up
Rscript Session-Info.R --pkglist-path $ppath
rm $ppath
##Run job stats script if it exists
command -v crc-job-stats &> /dev/null && command crc-job-stats
