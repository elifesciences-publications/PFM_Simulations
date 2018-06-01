#!/bin/sh

# Runs melodic and PROFUMO on simulated dataset
# Melodic is commited to queue, but PROFUMO is run directly (in nice), so make sure we're on Jalapeno18
# Usage: ICA_PROFUMO <TR> <filename> <dimensionality> <home directory>

tr=$1
Option=$2
dim=$3
homedir=$4

# Create file list for melodic
rm ${homedir}/Results/input_filelist.txt
for filename in `ls ${homedir}/Results/Temp/*_${Option}.nii.gz` ; do 
	echo $filename >> ${homedir}/Results/input_filelist.txt
done

# run melodic
fsl_sub -l logfiles -N M_${Option} melodic -i ${homedir}/Results/input_filelist.txt  \
-o ${homedir}/Results/Melodic_${Option}.gica --tr=${tr} --nobet --nomask -a concat --disableMigp -d ${dim}

# run PROFUMO (need to make sure we're on jalapeno18)
nice -n 20 ~samh/bin/PROFUMO ${homedir}/Results/PFMsims_atlas_${Option}.json ${dim} \
${homedir}/Results/PROFUMO_PFMsims_atlas_${Option} --useHRF ${tr} \
--hrfFile ~samh/PROFUMO/Scripts/DefaultHRF.phrf -d 0.05 > ${homedir}/Results/Output_PROFUMO_PFMsims_atlas_${Option}.txt

