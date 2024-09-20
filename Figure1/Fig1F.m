%% Detrend rpoD4 signal using basal expression level for each movie
clear 

load("WT_ML_24_DT_data.mat")

for i = 2:20
    xmins = []; 
    tmins = [];
    for j = 1:length(s_ML{i,2})
        x = s_ML{i,2}(j).MYs; 
        minx = min(x); 
        indx = find(x == minx); 
        xmins = cat(2,xmins,minx);
        tmins = cat(2,tmins,s_ML{i,2}(j).hrst(indx));
    end

    tmax = max([s_ML{i,2}.hrst]);

    xmins(tmins > tmax - 6) = [];
    tmins(tmins > tmax - 6) = [];

f = polyfit(tmins,xmins,2);

    for j = 1:length(s_ML{i,2})
        ts = s_ML{i,2}(j).hrst; 
        ys = s_ML{i,2}(j).MYs; 
        s_ML{i,2}(j).MYs_det = ys - ts.^2*f(1) - ts*f(2) - f(3);
    end
end
%%
s = s_ML;

pkmatun_wt = [];

for i = 2:20
    [pkmatun] = findpeaksmethod_rpoD4(s{i,2},'hrst','MYs_det','lengthMicrons_smooth');
    
    pkmatun_wt = cat(1,pkmatun_wt,pkmatun);
end



%% WT
%---- Peak amplitude of detrended data -----%
Amp = pkmatun_wt(:,5); % Peak amplitude of detrended data
CT = pkmatun_wt(:,4); % Circadian time

clear xs1 yperx1 lperx1

xs1 = unique(CT);
yperx1 = cell(length(xs1),1);
lperx1 = cell(length(xs1),1);

 for j = 1:length(CT)
        whereinxs = find(xs1 == CT(j));
        yperx1{whereinxs} = cat(1,yperx1{whereinxs},Amp(j));
 end

clear meanyperx1 stdyperx1 serrmean1

for i = 1:length(yperx1)
    meanyperx1(i) = mean(yperx1{i});
    stdyperx1(i) = std(yperx1{i});
    serrmean1(i) = stdyperx1(i)./sqrt(length(yperx1{i}));
end

clear xs2 meanyperx2 stdyperx2

xs2 = cat(1,xs1,xs1+24.75);
meanyperx2 = [meanyperx1,meanyperx1];
stdyperx2 = [stdyperx1,stdyperx1];

figure
curve1 = meanyperx2 + stdyperx2;
curve2 = meanyperx2 - stdyperx2;
x2 = [xs2', fliplr(xs2')];
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween,'c','LineStyle','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',0.2)
hold on;
% plot(xs2, meanyperx2, '-','color',[0.8500, 0.3250, 0.0980],'LineWidth', 2);
plot(xs2, meanyperx2, '-','color','c','LineWidth', 2);
%hold on;plot(xs2,smooth(meanyperx2,9),'b-','linewidth',2)
set(gca,'XTick',12*(0:40));
grid on
nightmodeon(s_ML{2,2},12,1);
xlim([0 48])
ylim([100 450])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Pulse amplitude (a.u.)','FontWeight','bold')
title('WT')
set(gca,'FontSize',18)
