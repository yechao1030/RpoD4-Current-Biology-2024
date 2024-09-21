%% Violin plot of birth length (Fig4C)
clear

load("Ptrc-rpoD4-2021-04-26_data.mat")

figure;
LB = [LB_WT;LB_rpoD4OX];
groupIdx = [ones(size(LB_WT));2*ones(size(LB_rpoD4OX))];
violin = violinplot(LB,groupIdx);
xticklabels({'WT','\itrpoD4OX'});
ylabel('Birth length (\mum)','fontweight','bold')
%ylim([2 5.5])
set(gca,'FontSize',16);

violin(1,1).ViolinColor = 'c';
violin(1,2).ViolinColor = 'm';
violin(1,1).ScatterPlot.SizeData = 10;
violin(1,2).ScatterPlot.SizeData = 10;


%% Violin plot of doubling time (Fig.S4C)
figure;
AGE = [AGE_WT;AGE_rpoD4OX];
groupIdx = [ones(size(AGE_WT));2*ones(size(AGE_rpoD4OX))];
violin = violinplot(AGE,groupIdx);
xticklabels({'WT','\itrpoD4OX'});
ylabel('Cell cycle duration (h)','fontweight','bold')
set(gca,'FontSize',16);

violin(1,1).ViolinColor = 'c';
violin(1,2).ViolinColor = 'm';
violin(1,1).ScatterPlot.SizeData = 10;
violin(1,2).ScatterPlot.SizeData = 10;


%% Violin plot of elongation rate (Fig.S4D)
figure;
ER = [ER_WT;ER_rpoD4OX];
groupIdx = [ones(size(ER_WT));2*ones(size(ER_rpoD4OX))];
violin = violinplot(ER,groupIdx);
xticklabels({'WT','\itrpoD4OX'});
ylabel('Elongation rate (h^{-1})','fontweight','bold')
set(gca,'FontSize',16);

violin(1,1).ViolinColor = 'c';
violin(1,2).ViolinColor = 'm';
violin(1,1).ScatterPlot.SizeData = 10;
violin(1,2).ScatterPlot.SizeData = 10;