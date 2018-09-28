#!/bin/sh

tr=0.72
dim=20
homedir='/home/fs0/janineb/scratch/HCP/DMN1200/PFM_Simulations/'
Option='MedO_MedM_MedT_01'

nice -n 20 ~samh/PROFUMO/PROFUMO_strong ${homedir}/Results/PFMsims_atlas_${Option}.json ${dim} \
${homedir}/Results/PROFUMO_PFMsims_atlas_${Option}_strong20 --useHRF ${tr} \
--hrfFile ~samh/PROFUMO/Scripts/DefaultHRF.phrf -d 0.05 > ${homedir}/Results/Output_PROFUMO_PFMsims_atlas_${Option}_strong20.txt

nice -n 20 ~samh/PROFUMO/PROFUMO_med ${homedir}/Results/PFMsims_atlas_${Option}.json ${dim} \
${homedir}/Results/PROFUMO_PFMsims_atlas_${Option}_med20 --useHRF ${tr} \
--hrfFile ~samh/PROFUMO/Scripts/DefaultHRF.phrf -d 0.05 > ${homedir}/Results/Output_PROFUMO_PFMsims_atlas_${Option}_med20.txt

nice -n 20 ~samh/PROFUMO/PROFUMO_weak ${homedir}/Results/PFMsims_atlas_${Option}.json ${dim} \
${homedir}/Results/PROFUMO_PFMsims_atlas_${Option}_weak20 --useHRF ${tr} \
--hrfFile ~samh/PROFUMO/Scripts/DefaultHRF.phrf -d 0.05 > ${homedir}/Results/Output_PROFUMO_PFMsims_atlas_${Option}_weak20.txt
