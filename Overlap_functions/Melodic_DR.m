function [sicaPg_new,sica_P1_DR_new,sica_A1_DR_new] = Melodic_DR(filename,D,atlasParams,params)

% Load existing data and melodic group map
sicaPg_new = read_avw(sprintf('Results/Melodic_%s.gica/melodic_IC.nii.gz',filename));
sicaPg_new = reshape(sicaPg_new,atlasParams.V,params.iN);

% Run stage 1 DR
sica_P1_DR_new = cell(params.S,1);
sica_A1_DR_new = cell(params.S,1); 
addpath ~steve/NETWORKS/FSLNets;
pinv_sicaPg_new = pinv(nets_demean(sicaPg_new)); 
for s = 1:params.S
    sica_A1_DR_new{s} = cell(params.R(s),1);
    M = zeros(atlasParams.V,params.iN,params.R(1));
    for r = 1:params.R(1)
        sica_A1_DR_new{s}{r} = nets_demean((pinv_sicaPg_new * nets_demean(double(D{s}{r})))');
        M(:,:,r) = demean((pinv(demean(double(sica_A1_DR_new{s}{r})))*demean(double(D{s}{r}))')');
    end
    sica_P1_DR_new{s} = mean(M,3);
end
clear s M r pinv_sicaPg_new






