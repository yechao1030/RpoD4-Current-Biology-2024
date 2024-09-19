%% Detrend rpoD4 signal using basal expression level for each movie

load("Fig3D_data.mat")

for i = 2:size(s,1)
    xmins = []; 
    tmins = [];
    for j = 1:length(s{i,2})
        x = s{i,2}(j).MYs; 
        minx = min(x); 
        indx = find(x == minx); 
        xmins = cat(2,xmins,minx);
        tmins = cat(2,tmins,s{i,2}(j).hrst(indx));
    end

    tmax = max([s{i,2}.hrst]);

    xmins(tmins > tmax - 6) = [];
    tmins(tmins > tmax - 6) = [];

f = polyfit(tmins,xmins,2);

    
    for j = 1:length(s{i,2})
        ts = s{i,2}(j).hrst; 
        ys = s{i,2}(j).MYs; 
        s{i,2}(j).MYs_det = ys - ts.^2*f(1) - ts*f(2) - f(3);
    end
end
%% Find rpoD4 expression peaks in the KaiCR393C mutant

pkmatun_wt = [];

for i = 2:size(s,1)
    [pkmatun] = findpeaksmethod_rpoD4_R393C(s{i,2},'hrst','MYs_det','lengthMicrons_smooth');
    
    pkmatun_wt = cat(1,pkmatun_wt,pkmatun);
end


%% KaiBCR393C
%---- Peak amplitude of detrended data normalised by mean pulse height -----%
Amp = pkmatun_wt(:,5); % Peak amplitude of detrended data
CT = pkmatun_wt(:,4); % Circadian time

clear xs1 yperx1 lperx1

xs1 = unique(CT);
yperx1 = cell(length(xs1),1);


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

meany = mean(meanyperx1);

clear xs2 meanyperx2 stdyperx2

xs2 = cat(1,xs1,xs1+16.5,xs1+16.5*2);
meanyperx2 = [meanyperx1,meanyperx1,meanyperx1];
stdyperx2 = [stdyperx1,stdyperx1,stdyperx1];

clear curve1 curve2 inBetween h 

figure;
curve1 = (meanyperx2 + stdyperx2)./meany;
curve2 = (meanyperx2 - stdyperx2)./meany;
x2 = [xs2', fliplr(xs2')];
inBetween = [curve1, fliplr(curve2)];
%h = fill(x2, inBetween,[0.8500, 0.3250, 0.0980],'LineStyle','none');
h = fill(x2, inBetween,[0.85 0.33 0.1],'LineStyle','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',0.2)
hold on;
% plot(xs2, meanyperx2, '-','color',[0.8500, 0.3250, 0.0980],'LineWidth', 2);
plot(xs2, meanyperx2./meany, '-','color',[0.85 0.33 0.1],'LineWidth', 2);
set(gca,'XTick',7.5*(0:40));
grid on
nightmodeon15(s{22,2},7.5,1);
xlim([0 30])
ylim([0.5 2.5])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Norm. pulse amplitude','FontWeight','bold')
%title('WT')
set(gca,'FontSize',18)



%% Plot Fig3D (no bining normalised)

figure;
plot(xs2, meanyperx2/meany, '-','color',[0.8500 0.3250 0.0980],'LineWidth', 2);
nightmodeon15(s{22,2},7.5,1);
xlim([0 30])
ylim([0.5 2])
set(gca,'XTick',7.5*(0:40));
grid on
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Norm. pulse amplitude','FontWeight','bold')
set(gca,'FontSize',14)





