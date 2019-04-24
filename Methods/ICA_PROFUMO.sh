#!/bin/sh

# Runs melodic and PROFUMO on simulated dataset
# Melodic is commited to queue, but PROFUMO is run directly (in nice), so make sure we're on Jalapeno18
# Usage: ICA_PROFUMO <basefilename> <nifti_dir> <dim> <TR>

basefilename=$1
nifti_dir=$2
dim=$3
TR=$4

# Create file list for melodic
for filename in ${nifti_dir}/*.nii.gz ; do
    echo $filename >> ${nifti_dir}/MELODIC_SpecFile.txt
done

# run melodic
fsl_sub -l logfiles -N M_$(basename $basefilename) \
    melodic -i ${nifti_dir}/MELODIC_SpecFile.txt \
    -o ${basefilename}_MELODIC.gica \
    --tr=${TR} --nobet --nomask -a concat -d ${dim}

# run PROFUMO (need to make sure we're on jalapeno18)
nice -n 20 ~samh/bin/PROFUMO \
    ${nifti_dir}/PROFUMO_SpecFile.json \
    ${dim} \
    ${basefilename}_PROFUMO.pfm \
    --useHRF ${TR} --hrfFile ~samh/PROFUMO/Scripts/DefaultHRF.phrf \
    --nThreads 15 -d 0.5 > ${basefilename}_PROFUMO_Output.txt

