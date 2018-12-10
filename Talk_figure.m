clear all; close all; clc

load('Results/Paper_figure_results.mat');

figure
histogram(Tnets_subs(:,2,:),'facecolor','g','edgecolor','g');
legend({'ICA'},'Location','NorthWest')
title('Accuracy of functional connectivity estimates')
xlabel(sprintf('Bad                                          Good'))
axis([ -0.2 1 0 180])
set(gca,'FontSize',16,'xtick',0:1)
print(gcf,'-dpng','-r300','Results/NAF2_talk1.png')

hold on;
histogram(Tnets_subs(:,1,:),'facecolor','r','edgecolor','r');
legend({'ICA','PFM'},'Location','NorthWest')
print(gcf,'-dpng','-r300','Results/NAF2_talk2.png')