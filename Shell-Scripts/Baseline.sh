#!/usr/bin/env bash

#SBATCH --job-name=Baseline
#SBATCH --qos=short
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=3
#SBATCH --output=../Outputs/Shell-Scripts/Baseline.out

bash Run-UMS.sh 0.0

