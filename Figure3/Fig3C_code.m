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
