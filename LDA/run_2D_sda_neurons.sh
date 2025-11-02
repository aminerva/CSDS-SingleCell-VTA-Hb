#!/bin/bash

#SBATCH --job-name=2D_sda_cv
#SBATCH --mail-user=aminerva@princeton.edu
#SBATCH --mail-type=START,END,FAIL
#SBATCH --time=3:00:00
#SBATCH --mem=240000
#SBATCH --cpus-per-task=4
#SBATCH --output='/jukebox/pena/Addie/scRNA-seq/CSDS/scripts/logs/%x_%j.log'

source /usr/people/aminerva/miniconda3/bin/activate /usr/people/aminerva/miniconda3/envs/scAnalysis 

# $1 = the cluster of interest

# Rscript --vanilla run_2D_sda_neurons.R "$1"
Rscript --vanilla run_2D_sda_neurons_cv_by_sample.R "$1"
