clear all; close all; clc

% Load & reshape results
load('Results/Paper_figure_results_New_sims_NoMisalignment.mat','Snets_subs','Tnets_subs');
Snets_NoM = shiftdim(Snets_subs,1); Tnets_NoM = shiftdim(Tnets_subs,1); clear Snet_subs Tnet_subs
load('Results/Paper_figure_results_New_sims_NoOverlap.mat','Snets_subs','Tnets_subs');
Snets_NoO = shiftdim(Snets_subs,1); Tnets_NoO = shiftdim(Tnets_subs,1); clear Snet_subs Tnet_subs
load('Results/Paper_figure_results_New_sims_HighM.mat','Snets_subs','Tnets_subs');
Snets = shiftdim(Snets_subs,1); Tnets = shiftdim(Tnets_subs,1); clear Snet_subs Tnet_subs
Snets_NoM = reshape(Snets_NoM,size(Snets_NoM,1),size(Snets_NoM,2)*size(Snets_NoM,3))';
Tnets_NoM = reshape(Tnets_NoM,size(Tnets_NoM,1),size(Tnets_NoM,2)*size(Tnets_NoM,3))';
Snets_NoO = reshape(Snets_NoO,size(Snets_NoO,1),size(Snets_NoO,2)*size(Snets_NoO,3))';
Tnets_NoO = reshape(Tnets_NoO,size(Tnets_NoO,1),size(Tnets_NoO,2)*size(Tnets_NoO,3))';
Snets = reshape(Snets,size(Snets,1),size(Snets,2)*size(Snets,3))';
Tnets = reshape(Tnets,size(Tnets,1),size(Tnets,2)*size(Tnets,3))';

% Plot
Nmethods = {'Traditional ICA-DR','Thresholded ICA-DR','PROFUMO'};
Ntests = {'full simulation','no misalignment','no overlap'};
color = ['g','b','r','g','b','r','g','b','r'];

figure; set(gcf,'Position',[0 0 650 550],'PaperPositionMode','auto')
subplot(2,1,1);
X = [Tnets(:,2); Tnets(:,3); Tnets(:,1); Tnets_NoM(:,2); Tnets_NoM(:,3); Tnets_NoM(:,1); Tnets_NoO(:,2); Tnets_NoO(:,3); Tnets_NoO(:,1)];
G1 = [ones(size(Tnets,1)*3,1); repmat(2,size(Tnets,1)*3,1); repmat(3,size(Tnets,1)*3,1)];
G2 = [ones(size(Tnets,1),1); ones(size(Tnets,1),1)+1; ones(size(Tnets,1),1)+2; ones(size(Tnets,1),1); ones(size(Tnets,1),1)+1; ones(size(Tnets,1),1)+2; ones(size(Tnets,1),1); ones(size(Tnets,1),1)+1; ones(size(Tnets,1),1)+2];
G1(isnan(X)==1) = []; G2(isnan(X)==1) = []; X(isnan(X)==1) = []; 
boxplot(X,{G1, G2},'factorgap',10,'colors',color,'symbol','ko','outliersize',1);
set(gca,'xtick',[2 5.8 9.6],'xticklabel',Ntests)
hold on; hline(0,':k')
title('Accuracy of temporal edge estimation'); 
subplot(2,1,2);
X = [Snets(:,2); Snets(:,3); Snets(:,1); Snets_NoM(:,2); Snets_NoM(:,3); Snets_NoM(:,1); Snets_NoO(:,2); Snets_NoO(:,3); Snets_NoO(:,1)];
G1 = [ones(size(Tnets,1)*3,1); repmat(2,size(Tnets,1)*3,1); repmat(3,size(Tnets,1)*3,1)];
G2 = [ones(size(Tnets,1),1); ones(size(Tnets,1),1)+1; ones(size(Tnets,1),1)+2; ones(size(Tnets,1),1); ones(size(Tnets,1),1)+1; ones(size(Tnets,1),1)+2; ones(size(Tnets,1),1); ones(size(Tnets,1),1)+1; ones(size(Tnets,1),1)+2];
G1(isnan(X)==1) = []; G2(isnan(X)==1) = []; X(isnan(X)==1) = []; 
boxplot(X,{G1, G2},'factorgap',10,'colors',color,'symbol','ko','outliersize',1);
set(gca,'xtick',[2 5.8 9.6],'xticklabel',Ntests)
hold on; hline(0,':k')
title('Accuracy of spatial edge estimation');
print(gcf,'-dpng','-r300','Results/Paper_figure_comparisons.png')

% Do stats
Snets = 0.5*log((1+Snets)./(1-Snets)); Tnets = 0.5*log((1+Tnets)./(1-Tnets));
Snets_NoM = 0.5*log((1+Snets_NoM)./(1-Snets_NoM)); Tnets_NoM = 0.5*log((1+Tnets_NoM)./(1-Tnets_NoM));
Snets_NoO = 0.5*log((1+Snets_NoO)./(1-Snets_NoO)); Tnets_NoO = 0.5*log((1+Tnets_NoO)./(1-Tnets_NoO));

