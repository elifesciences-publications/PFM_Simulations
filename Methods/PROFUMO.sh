#!/bin/sh

# Runs PROFUMO on simulated dataset
# Usage: PROFUMO.sh <basefilename> <nifti_dir> <dim> <TR>

basefilename=$1
nifti_dir=$2
dim=$3
TR=$4

# Run PROFUMO (need to make sure we're on jalapeno18!)
nice -n 20 ~samh/bin/PROFUMO \
    ${nifti_dir}/PROFUMO_SpecFile.json \
    ${dim} \
    ${basefilename}_PROFUMO.pfm \
    --useHRF ${TR} --hrfFile ~samh/PROFUMO/Scripts/DefaultHRF.phrf \
    --nThreads 15 -d 0.5 \
    > ${basefilename}_PROFUMO_Output.txt
