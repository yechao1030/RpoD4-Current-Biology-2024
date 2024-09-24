%% Fig.S3A
clear

load("FigS3_data.mat")

WT_YFP_mean = nanmean(WT_YFP_s);
WT_YFP_std = nanstd(WT_YFP_s);

delD4_YFP_mean = nanmean(delD4_YFP_s);
delD4_YFP_std = nanstd(delD4_YFP_s);

xs = t1;

figure
plot(xs,WT_YFP_mean,'c-','linewidth',2)
hold on;
inBetween = [WT_YFP_mean+WT_YFP_std, fliplr(WT_YFP_mean-WT_YFP_std)];
x2 = [xs, fliplr(xs)];
h = fill(x2, inBetween,'c','LineStyle','none');set(h,'facealpha',0.2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

clear inBetween

plot(xs,delD4_YFP_mean,'b-','linewidth',2)
inBetween = [delD4_YFP_mean+delD4_YFP_std, fliplr(delD4_YFP_mean-delD4_YFP_std)];
x2 = [xs, fliplr(xs)];
h = fill(x2, inBetween,'b','LineStyle','none');set(h,'facealpha',0.2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

nightmodeon(s_JRC35{2,2},12,1)
xlim([0 144])
ylim([300 800])
set(gca,'XTick',12*(0:40));
grid on;
xlabel('Time in LL (h)','fontweight','bold')
ylabel('P_{\itkaiBC}-EYFP (a.u.)','fontweight','bold')
legend({'WT','\it\DeltarpoD4'},'Location','northwest')

%% Fig.S3B
figure; % box plot of period (BioDare2, detrended data)
period = [WT_JRC35_dtr_period;delD4_JRC35_dtr_period];
groupIdx = [ones(size(WT_JRC35_dtr_period));2*ones(size(delD4_JRC35_dtr_period))];
plot_box_scatter(period,groupIdx,'color',{'c','b'},'symbol',{'o','o'},'opt',1);
xticklabels({'WT','ΔrpoD4'});
ylabel('Period (h)','fontweight','bold')
ylim([23.5 26.5])
set(gca,'FontSize',14);

%% Fig.S3C
rpoD4_YFP_mean = nanmean(rpoD4_YFP_s);
rpoD4_YFP_std = nanstd(rpoD4_YFP_s);

rpoD4OX_YFP_mean = nanmean(rpoD4OX_YFP_s);
rpoD4OX_YFP_std = nanstd(rpoD4OX_YFP_s);

xs = t2;

figure
plot(xs,rpoD4_YFP_mean,'k-','linewidth',2)
hold on;
inBetween = [rpoD4_YFP_mean+rpoD4_YFP_std, fliplr(rpoD4_YFP_mean-rpoD4_YFP_std)];
x2 = [xs, fliplr(xs)];
h = fill(x2, inBetween,'k','LineStyle','none');set(h,'facealpha',0.2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

clear inBetween

plot(xs,rpoD4OX_YFP_mean,'m-','linewidth',2)
inBetween = [rpoD4OX_YFP_mean+rpoD4OX_YFP_std, fliplr(rpoD4OX_YFP_mean-rpoD4OX_YFP_std)];
x2 = [xs, fliplr(xs)];
h = fill(x2, inBetween,'m','LineStyle','none');set(h,'facealpha',0.2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

nightmodeon(s_PtrcD4OX{6,2},12,1)
xlim([0 144])
ylim([500 1400])
set(gca,'XTick',12*(0:40));
grid on;
xlabel('Time in LL (h)','fontweight','bold')
ylabel('P_{\itkaiBC}-EYFP (a.u.)','fontweight','bold')
legend({'\itrpoD4 non-induced','\itrpoD4OX'},'Location','northwest')

%% Fig.S3D
figure; % box plot of period (BioDare2, detrended data)
period = [rpoD4_kaiBC_dtr_period;rpoD4OX_kaiBC_dtr_period];
groupIdx = [ones(size(rpoD4_kaiBC_dtr_period));2*ones(size(rpoD4OX_kaiBC_dtr_period))];
plot_box_scatter(period,groupIdx,'color',{'k','m'},'symbol',{'o','o'},'opt',1);
xticklabels({'rpoD4 non-induced','rpoD4OX'});
ylabel('Period (h)','fontweight','bold')
ylim([24 26.5])
set(gca,'FontSize',14);