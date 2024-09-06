%% Plot Figure 3C

load("Fig3C_data.mat"); % Here I load the proccessed data for plotting Fig3C directly.
% To generated these "processed" data, please read the instruction document
% in this folder. CY 2024.09.06

figure;
plot(xs2_VLL,meanyperx2_VLL/meany_VLL,'-','Color',[0.00,0.45,0.74],'LineWidth',2);
hold on;
plot(xs2_ML,meanyperx2_ML/meany_ML,'-','Color',[0.47,0.67,0.19],'LineWidth',2);
plot(xs2_HL,meanyperx2_HL/meany_HL,'-','Color',[0.85,0.33,0.10],'LineWidth',2);
set(gca,'XTick',12*(0:40));grid on
nightmodeon(s_VLL{5,2},12,1);
xlim([0 48])
ylim([0.5 2.5])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Norm. pulse amplitude','FontWeight','bold')
leg = legend({'10 \mumol m^{-2} s^{-1}','18 \mumol m^{-2} s^{-1}','25 \mumol m^{-2} s^{-1}'},'Location','northeast');
leg.ItemTokenSize = [20,5];

%%
figure;
plot(xs2_ML,meanyperx2_ML,'-','Color',[0.47,0.67,0.19],'LineWidth',2);
hold on
plot(xs2_HL,meanyperx2_HL,'-','Color',[0.85,0.33,0.10],'LineWidth',2);
set(gca,'XTick',12*(0:40));grid on
nightmodeon(s_ML{3,2},12,1);
xlim([0 48])
ylim([100 500])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Pulse amplitude (a.u.)','FontWeight','bold')
leg = legend({'18 \mumol m^{-2} s^{-1}','25 \mumol m^{-2} s^{-1}'},'Location','northeast');
leg.ItemTokenSize = [20,5];



%%
figure;
plot(xs2_VLL,meanyperx2_VLL./meany_VLL,'-','Color',[0.00,0.45,0.74],'LineWidth',2);
hold on
plot(xs2_ML,meanyperx2_ML./meany_ML,'-','Color',[0.47,0.67,0.19],'LineWidth',2);
hold on
plot(xs2_HL,meanyperx2_HL./meany_HL,'-','Color',[0.85,0.33,0.10],'LineWidth',2);
set(gca,'XTick',12*(0:40));grid on
nightmodeon(s_ML{3,2},12,1);
xlim([0 48])
ylim([0 2.5])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Norm. pulse amplitude','FontWeight','bold')
leg = legend({'10 \mumol m^{-2} s^{-1}','18 \mumol m^{-2} s^{-1}','25 \mumol m^{-2} s^{-1}'},'Location','northeast');
leg.ItemTokenSize = [20,5];

%% Smooth
figure;
plot(xs2_VLL./24.75*24,smooth(meanyperx2_VLL./meany_VLL,5),'-','Color',[0.00,0.45,0.74],'LineWidth',2);
hold on
plot(xs2_ML./24.75*24,smooth(meanyperx2_ML./meany_ML,5),'-','Color',[0.47,0.67,0.19],'LineWidth',2);
hold on
plot(xs2_HL./24.75*24,smooth(meanyperx2_HL./meany_HL,5),'-','Color',[0.85,0.33,0.10],'LineWidth',2);
set(gca,'XTick',12*(0:40));grid on
xlim([0 49.5])
ylim([0.5 2])
xlabel('Circadian time of day (h)','FontWeight','bold')
ylabel('Norm. pulse amplitude','FontWeight','bold')
leg = legend({'10 \mumol m^{-2} s^{-1}','18 \mumol m^{-2} s^{-1}','25 \mumol m^{-2} s^{-1}'},'Location','northeast');
leg.ItemTokenSize = [20,5];


%% binning 1.5 h normalised
meanyperx1_VLL_bin = [];
xs1_VLL_bin = [];

for i = 1:size(meanyperx2_VLL,2)/4
    meanyperx1_VLL_bin(i) = (meanyperx2_VLL(i*2-1)+meanyperx2_VLL(i*2))/2;
    xs1_VLL_bin(i) = (xs2_VLL(i*2-1)+xs2_VLL(i*2))/2;
end

xs2_VLL_bin = [xs1_VLL_bin,xs1_VLL_bin+24];
meanyperx2_VLL_bin = [meanyperx1_VLL_bin,meanyperx1_VLL_bin];

meanyperx1_ML_bin = [];
xs1_ML_bin = [];

for i = 1:size(meanyperx2_ML,2)/4
    meanyperx1_ML_bin(i) = (meanyperx2_ML(i*2-1)+meanyperx2_ML(i*2))/2;
    xs1_ML_bin(i) = (xs2_ML(i*2-1)+xs2_ML(i*2))/2;
end

xs2_ML_bin = [xs1_ML_bin,xs1_ML_bin+24];
meanyperx2_ML_bin = [meanyperx1_ML_bin,meanyperx1_ML_bin];

meanyperx1_HL_bin = [];
xs1_HL_bin = [];

for i = 1:size(meanyperx2_HL,2)/4
    meanyperx1_HL_bin(i) = (meanyperx2_HL(i*2-1)+meanyperx2_HL(i*2))/2;
    xs1_HL_bin(i) = (xs2_HL(i*2-1)+xs2_HL(i*2))/2;
end

xs2_HL_bin = [xs1_HL_bin,xs1_HL_bin+24];
meanyperx2_HL_bin = [meanyperx1_HL_bin,meanyperx1_HL_bin];

figure;
plot(xs2_VLL_bin,meanyperx2_VLL_bin./meany_VLL,'-','Color',[0.00,0.45,0.74],'LineWidth',2);
hold on
plot(xs2_ML_bin,meanyperx2_ML_bin./meany_ML,'-','Color',[0.47,0.67,0.19],'LineWidth',2);
hold on
plot(xs2_HL_bin,meanyperx2_HL_bin./meany_HL,'-','Color',[0.85,0.33,0.10],'LineWidth',2);
set(gca,'XTick',12*(0:40));grid on
nightmodeon(s_ML{3,2},12,1);
xlim([0 48])
ylim([0.5 2.5])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Norm. pulse amplitude','FontWeight','bold')
leg = legend({'10 \mumol m^{-2} s^{-1}','18 \mumol m^{-2} s^{-1}','25 \mumol m^{-2} s^{-1}'},'Location','northeast');
leg.ItemTokenSize = [20,5];