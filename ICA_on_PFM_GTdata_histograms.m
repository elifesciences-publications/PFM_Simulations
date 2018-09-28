clear all; close all; clc

load('temp.mat')
for n = 2:2:size(ICAall,2)
    figure; set(gcf,'Position',[0 0 1800 900],'PaperPositionMode','auto')
    subplot(2,1,1);
    histogram(ICAall(:,n-1),-5:0.05:5); hold on; histogram(PFMall(:,n-1),-5:0.05:5);
    subplot(2,1,2);
    histogram(ICAall(:,n),-5:0.05:5); hold on; histogram(PFMall(:,n),-5:0.05:5);
    legend({'ICA','PFM'})
    print(sprintf('%02d_example_hist',n/2),'-dpng','-r150'); 
end
   
clear all; close all; clc

load('temp_real.mat')
for n = 2:2:size(ICAall,2)
    figure; set(gcf,'Position',[0 0 1800 900],'PaperPositionMode','auto')
    subplot(2,1,1);
    histogram(ICAall(:,n-1),-5:0.05:5); hold on; histogram(PFMall(:,n-1),-5:0.05:5);
    subplot(2,1,2);
    histogram(ICAall(:,n),-5:0.05:5); hold on; histogram(PFMall(:,n),-5:0.05:5);
    legend({'ICA','PFM'})
    print(sprintf('%02d_example_real_hist',n/2),'-dpng','-r150'); 
end

