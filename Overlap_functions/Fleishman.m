function Y = Fleishman(X)

% Load real data skew & kurtosis (based on 100307/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii)
% Also load Fleishman power method coefficient table ('https://chandrakant721.wordpress.com/fleishmans-power-method-coefficient-table/')
load('Overlap_functions/Fleish.mat');

% Pick random kurtosis & skewness pair
I = randi(size(sk,1),1);
Iskew = round(abs(sk(I)),1);
Ikurt = round(ku(I),1);
sk_match = find(F(:,1)==Iskew);
ku_match = find(F(:,2)==Ikurt);
all_match = intersect(sk_match,ku_match);
b = F(all_match,3);
c = F(all_match,4);
d = F(all_match,5);

% Use third polynomial approach to add the non-Gaussianity
% Based on ('https://chandrakant721.wordpress.com/2015/10/25/simulation-of-non-normal-distribution-in-simple-steps/')
Y = -c + b*X + c*X.^2 + d*X.^3; 

