#!/bin/sh

# Runs melodic and PROFUMO on simulated dataset
# Melodic is commited to queue, but PROFUMO is run directly (in nice), so make sure we're on Jalapeno18
# Usage: ICA_PROFUMO <TR> <filename> <dimensionality> <home directory>
# Usage: ICA_PROFUMO <filename> <nifti_dir> <dim> <TR>

filename=$1
nifti_dir=$2
dim=$3
TR=$4

# Create file list for melodic
for filename in ${nifti_dir}/*.nii.gz ; do
    echo $filename >> ${nifti_dir}/MELODIC_SpecFile.txt
done

# run melodic
fsl_sub -l logfiles -N M_${Option} \
    melodic -i ${nifti_dir}/MELODIC_SpecFile.txt \
    -o ${filename}_MELODIC.gica \
    --tr=${TR} --nobet --nomask -a concat --disableMigp -d ${dim}

# run PROFUMO (need to make sure we're on jalapeno18)
nice -n 20 ~samh/bin/PROFUMO \
    ${nifti_dir}/PROFUMO_SpecFile.txt \
    ${dim} \
    ${filename}_PROFUMO.pfm \
    --useHRF ${TR} --hrfFile ~samh/PROFUMO/Scripts/DefaultHRF.phrf \
    --nThreads 15 -d 0.5 > ${filename}_PROFUMO_Output.txt

