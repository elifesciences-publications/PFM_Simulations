#!/bin/sh

# Runs MELODIC on simulated dataset
# Usage: MELODIC.sh <basefilename> <nifti_dir> <dim> <TR>

basefilename=$1
nifti_dir=$2
dim=$3
TR=$4

# Run MELODIC
# fsl_sub -l logfiles -N M_$(basename $basefilename) \
melodic -i ${nifti_dir}/MELODIC_SpecFile.txt \
    -o ${basefilename}_MELODIC.gica \
    --tr=${TR} -a concat -d ${dim} \
    --migpN=300 --nobet --nomask \
    --maxrestart=5 --verbose \
    > ${basefilename}_MELODIC_Output.txt
