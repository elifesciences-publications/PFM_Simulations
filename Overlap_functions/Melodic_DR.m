function [sicaPg_new,sica_P1_DR_new,sica_A1_DR_new,sica_P1_DR_thres,sica_A1_DR_thres_noO,sica_A1_DR_thres_withO] = Melodic_DR(filename,D,atlasParams,params,aggressive)

addpath ~steve/NETWORKS/FSLNets;
addpath ~/scratch/matlab/;

% Which thresholds to use
if aggressive == 1
    THpos = 1; THneg = 4;
elseif aggressive == 0
    THpos = 3; THneg = 2;
end

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

% Run triple regression after thresholding
sica_P1_DR_thres = cell(params.S,1);
sica_A1_DR_thres_noO = cell(params.S,1);
sica_A1_DR_thres_withO = cell(params.S,1);

for s = 1:params.S
    M = sica_P1_DR_new{s};
        
    % Apply mixture modelling to each map to threshold:
    M1 = zeros(size(M));
    for n = 1:size(M,2)
        [~,thresh] = ggfit(M(:,n));
        M1(M(:,n)<thresh(THneg),n) = M(M(:,n)<thresh(THneg),n);
        M1(M(:,n)>thresh(THpos),n) = M(M(:,n)>thresh(THpos),n);
    end
    sica_P1_DR_thres{s} = M1;
    
    % Remove overlap (again using mixture modelling to threshold):
    M = M1;
    O = sum(abs(M1),2); 
    [~,thresh] = ggfit(O);
    M1(O>thresh(3),:) = 0;
    
    % If a map is missing entirely, add it back in from the thresholded version:
    if find(sum(M1)==0); M1(:,sum(M1)==0) = M(:,sum(M1)==0); end
    
    % Obtain timeseries via simple masking:
    sica_A1_DR_thres_noO{s} = cell(params.R(s),1);
    sica_A1_DR_thres_withO{s} = cell(params.R(s),1);
    for r = 1:params.R(1)
        %sica_A1_DR_thres{s}{r} = nets_demean((pinv(nets_demean(M)) * nets_demean(double(D{s}{r})))');
        sica_A1_DR_thres_noO{s}{r} = nets_demean((M1'*nets_demean(double(D{s}{r})))');
        sica_A1_DR_thres_withO{s}{r} = nets_demean((pinv(nets_demean(M)) * nets_demean(double(D{s}{r})))');
    end
end

% Clear variables
clear s M r 

%%%%%%%%%%%%%% OLD OVERLAP STUFF %%%%%%%%%%%%%%%

%     % normalise and fix sign
%     M1 = M.*repmat(sign(mean(M)),size(M,1),1)./repmat(max(abs(M)),size(M,1),1);
%     [C12,~] = spatialcorr(M1,M); C12 = sign(C12(eye(size(C12,1))==1));
%     M1 = M1.*repmat(C12',size(M1,1),1); 
%     % threshold
%     M1(M1<0.3) = 0;
%     % remove overlap
%     for n = 1:params.N
%         if n == 1; O = M1(:,1); else O = O.*M1(:,n); end
%     end
%     M1(O>0.2,:) = 0;
%     %figure; subplot(1,2,1); imagesc(M,[-40 40]); colorbar; subplot(1,2,2); imagesc(M1,[-1 1]); colorbar  







