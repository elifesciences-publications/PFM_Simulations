function [ Pg ] = loadMELODIC( baseFileName, params )

Pg = read_avw([baseFileName '_MELODIC.gica/melodic_IC.nii.gz']);
Pg = reshape(Pg, params.V, params.iN);

return
