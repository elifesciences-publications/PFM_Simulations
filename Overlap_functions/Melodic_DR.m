function [sicaPg_new,sica_P1_DR_new,sica_A1_DR_new,sica_P1_DR_thres,sica_A1_DR_thres] = Melodic_DR(filename,D,atlasParams,params)

addpath ~steve/NETWORKS/FSLNets;

% Load existing data and melodic group map
sicaPg_new = read_avw(sprintf('Results/Melodic_%s.gica/melodic_IC.nii.gz',filename));
sicaPg_new = reshape(sicaPg_new,atlasParams.V,params.iN);

% Run DR
sica_P1_DR_new = cell(params.S,1);
sica_A1_DR_new = cell(params.S,1); 
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

% Run DR after thresholding
sica_P1_DR_thres = cell(params.S,1);
sica_A1_DR_thres = cell(params.S,1); 
sicaPg_new_thres = sicaPg_new;
sicaPg_new_thres(abs(sicaPg_new_thres)<5) = 0;
pinv_sicaPg_new = pinv(nets_demean(sicaPg_new_thres)); 
for s = 1:params.S
    sica_A1_DR_thres{s} = cell(params.R(s),1);
    M = zeros(atlasParams.V,params.iN,params.R(1));
    for r = 1:params.R(1)
        sica_A1_DR_thres{s}{r} = nets_demean((pinv_sicaPg_new * nets_demean(double(D{s}{r})))');
        M(:,:,r) = demean((pinv(demean(double(sica_A1_DR_thres{s}{r})))*demean(double(D{s}{r}))')');
    end
    sica_P1_DR_thres{s} = mean(M,3);
end

% Clear variables
clear s M r pinv_sicaPg_new






