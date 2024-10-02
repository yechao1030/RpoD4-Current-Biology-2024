%% Get unique rpoD4 expression peak information from time traces (10 umol)
s = s_VLL;
pkmatun_wt = [];

for i = 2:11
    [pkmatun] = findpeaksmethod_rpoD4(s{i,2},'hrst','MYs_det','lengthMicrons_smooth',24,0,120,25);
    pkmatun_wt = cat(1,pkmatun_wt,pkmatun);
end


%% WT
%---- Peak amplitude of detrended data normalised by mean pulse height -----%
Amp = pkmatun_wt(:,5); % Peak amplitude of detrended data
CT = pkmatun_wt(:,4); % Circadian time

clear xs1 yperx1 xs2 meanyperx1 stdyperx1 meanyperx2 stdyperx2 s2

xs1 = unique(CT);
yperx1 = cell(length(xs1),1);


 for j = 1:length(CT)
        whereinxs = find(xs1 == CT(j));
        yperx1{whereinxs} = cat(1,yperx1{whereinxs},Amp(j));
 end


for i = 1:length(yperx1)
    meanyperx1(i) = mean(yperx1{i});
    stdyperx1(i) = std(yperx1{i});
    serrmean1(i) = stdyperx1(i)./sqrt(length(yperx1{i}));
end

meany_VLL = mean(meanyperx1);

xs2 = cat(1,xs1,xs1+24);%24-26
meanyperx2 = [meanyperx1,meanyperx1];
stdyperx2 = [stdyperx1,stdyperx1];

figure;
curve1 = (meanyperx2 + stdyperx2)./meany_VLL;
curve2 = (meanyperx2 - stdyperx2)./meany_VLL;
x2 = [xs2', fliplr(xs2')];
inBetween = [curve1, fliplr(curve2)];
%h = fill(x2, inBetween,[0.8500, 0.3250, 0.0980],'LineStyle','none');
h = fill(x2, inBetween,'c','LineStyle','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',0.2)
hold on;
% plot(xs2, meanyperx2, '-','color',[0.8500, 0.3250, 0.0980],'LineWidth', 2);
plot(xs2, meanyperx2./meany_VLL, '-','color','c','LineWidth', 2);
set(gca,'XTick',12*(0:40));
grid on
nightmodeon(s_VLL{10,2},12,1);
xlim([0 48])
ylim([0 2.5])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Norm. pulse height (a.u.)','FontWeight','bold')
%title('WT')
set(gca,'FontSize',18)

meanyperx2_VLL = meanyperx2;
xs2_VLL = xs2;


%% WT
%---- Peak amplitude of non-detrended data -----%
Amp = pkmatun_wt(:,5); % Peak amplitude of detrended data
CT = pkmatun_wt(:,4); % Circadian time

clear xs1 yperx1 xs2 meanyperx1 stdyperx1 meanyperx2 stdyperx2 s2

xs1 = unique(CT);
yperx1 = cell(length(xs1),1);


 for j = 1:length(CT)
        whereinxs = find(xs1 == CT(j));
        yperx1{whereinxs} = cat(1,yperx1{whereinxs},Amp(j));
 end


for i = 1:length(yperx1)
    meanyperx1(i) = mean(yperx1{i});
    stdyperx1(i) = std(yperx1{i});
    serrmean1(i) = stdyperx1(i)./sqrt(length(yperx1{i}));
end

xs2 = cat(1,xs1,xs1+24);%24-26
meanyperx2 = [meanyperx1,meanyperx1];
stdyperx2 = [stdyperx1,stdyperx1];

figure;
curve1 = meanyperx2 + stdyperx2;
curve2 = meanyperx2 - stdyperx2;
x2 = [xs2', fliplr(xs2')];
inBetween = [curve1, fliplr(curve2)];
%h = fill(x2, inBetween,[0.8500, 0.3250, 0.0980],'LineStyle','none');
h = fill(x2, inBetween,'c','LineStyle','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',0.2)
hold on;
% plot(xs2, meanyperx2, '-','color',[0.8500, 0.3250, 0.0980],'LineWidth', 2);
plot(xs2, meanyperx2, '-','color','c','LineWidth', 2);
set(gca,'XTick',12*(0:40));
grid on
nightmodeon(s_VLL{10,2},12,1);
xlim([0 48])
ylim([100 250])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Pulse amplitude (a.u.)','FontWeight','bold')
%title('WT')
set(gca,'FontSize',18)


%%
%---- Peak width of non-detrended data -----%
Width = pkmatun_wt(:,6); % Peak width of detrended data
CT = pkmatun_wt(:,4); % Circadian time

clear xs1 yperx1 xs2 meanyperx1 stdyperx1 meanyperx2 stdyperx2 s2

xs1 = unique(CT);
yperx1 = cell(length(xs1),1);


 for j = 1:length(CT)
        whereinxs = find(xs1 == CT(j));
        yperx1{whereinxs} = cat(1,yperx1{whereinxs},Width(j));
 end

for i = 1:length(yperx1)
    meanyperx1(i) = mean(yperx1{i});
    stdyperx1(i) = std(yperx1{i});
    serrmean1(i) = stdyperx1(i)./sqrt(length(yperx1{i}));
end

xs2 = cat(1,xs1,xs1+24);
meanyperx2 = [meanyperx1,meanyperx1];
stdyperx2 = [stdyperx1,stdyperx1];

figure;
curve1 = meanyperx2 + stdyperx2;
curve2 = meanyperx2 - stdyperx2;
x2 = [xs2', fliplr(xs2')];
inBetween = [curve1, fliplr(curve2)];
%h = fill(x2, inBetween,[0.8500, 0.3250, 0.0980],'LineStyle','none');
h = fill(x2, inBetween,'c','LineStyle','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',0.2)
hold on;
% plot(xs2, meanyperx2, '-','color',[0.8500, 0.3250, 0.0980],'LineWidth', 2);
plot(xs2, meanyperx2, '-','color','c','LineWidth', 2);
set(gca,'XTick',12*(0:40));
grid on
nightmodeon(s_VLL{10,2},12,1);
xlim([0 48])
ylim([0 14])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Pulse width (h)','FontWeight','bold')
%title('WT')
set(gca,'FontSize',18)