%% WT (10 umol) non-detrended data --- start of Fig.S2A,B&C
clear

load ('FigS2_data.mat')

s = s_VLL;
pkmatun_wt = [];

for i = 2:11
    [pkmatun] = findpeaksmethod_rpoD4_v2(s{i,2},'hrst','MYs','lengthMicrons_smooth');
    pkmatun_wt = cat(1,pkmatun_wt,pkmatun);
end

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

xs2 = cat(1,xs1,xs1+24);
meanyperx2 = [meanyperx1,meanyperx1];
stdyperx2 = [stdyperx1,stdyperx1];

figure;
curve1 = meanyperx2 + stdyperx2;
curve2 = meanyperx2 - stdyperx2;
x2 = [xs2', fliplr(xs2')];
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween,'c','LineStyle','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',0.2)
hold on;
plot(xs2, meanyperx2, '-','color','c','LineWidth', 2);
set(gca,'XTick',12*(0:40));
grid on
nightmodeon(s_VLL{10,2},12,1);
xlim([0 48])
ylim([100 250])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Pulse amplitude (a.u.)','FontWeight','bold')
title('10 \mumol m^{-2} s^{-1}')
set(gca,'FontSize',18)

%%  WT (10 umol) detrended data
s = s_VLL;
pkmatun_wt = [];

for i = 2:11
    [pkmatun] = findpeaksmethod_rpoD4_v2(s{i,2},'hrst','MYs_det','lengthMicrons_smooth');
    pkmatun_wt = cat(1,pkmatun_wt,pkmatun);
end

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

xs2 = cat(1,xs1,xs1+24);
meanyperx2 = [meanyperx1,meanyperx1];
stdyperx2 = [stdyperx1,stdyperx1];

figure;
curve1 = meanyperx2 + stdyperx2;
curve2 = meanyperx2 - stdyperx2;
x2 = [xs2', fliplr(xs2')];
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween,'c','LineStyle','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',0.2)
hold on;
plot(xs2, meanyperx2, '-','color','c','LineWidth', 2);
set(gca,'XTick',12*(0:40));
grid on
nightmodeon(s_VLL{10,2},12,1);
xlim([0 48])
ylim([0 120])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Pulse amplitude (a.u.)','FontWeight','bold')
title('10 \mumol m^{-2} s^{-1}')
set(gca,'FontSize',18)

%% WT (10 umol) detrended and normalised data

s = s_VLL;
pkmatun_wt = [];

for i = 2:11
    [pkmatun] = findpeaksmethod_rpoD4_v2(s{i,2},'hrst','MYs_det','lengthMicrons_smooth');
    pkmatun_wt = cat(1,pkmatun_wt,pkmatun);
end

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
h = fill(x2, inBetween,'c','LineStyle','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',0.2)
hold on;
plot(xs2, meanyperx2./meany_VLL, '-','color','c','LineWidth', 2);
set(gca,'XTick',12*(0:40));
grid on
nightmodeon(s_VLL{10,2},12,1);
xlim([0 48])
ylim([0 2.5])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Norm. pulse height (a.u.)','FontWeight','bold')
title('10 \mumol m^{-2} s^{-1}')
set(gca,'FontSize',18)


%% WT (18 umol) non-detrended data
s = s_ML;

pkmatun_wt = [];

for i = 2:20
    [pkmatun] = findpeaksmethod_rpoD4(s{i,2},'hrst','MYs','lengthMicrons_smooth');
    
    pkmatun_wt = cat(1,pkmatun_wt,pkmatun);
end

%---- Peak amplitude of non-detrended data -----%
Amp = pkmatun_wt(:,5); % Peak amplitude of detrended data
CT = pkmatun_wt(:,4); % Circadian time

clear xs1 yperx1 xs2 meanyperx1 stdyperx1 meanyperx2 stdyperx2 s2

xs1 = unique(CT);
yperx1 = cell(length(xs1),1);
lperx1 = cell(length(xs1),1);

 for j = 1:length(CT)
        whereinxs = find(xs1 == CT(j));
        yperx1{whereinxs} = cat(1,yperx1{whereinxs},Amp(j));
 end

for i = 1:length(yperx1)
    meanyperx1(i) = mean(yperx1{i});
    stdyperx1(i) = std(yperx1{i});
    serrmean1(i) = stdyperx1(i)./sqrt(length(yperx1{i}));
end

xs2 = cat(1,xs1,xs1+24);
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
plot(xs2, meanyperx2, '-','color','c','LineWidth', 2);
set(gca,'XTick',12*(0:40));
grid on
nightmodeon(s_ML{2,2},12,1);
xlim([0 48])
ylim([600 1100])
title('18 \mumol m^{-2} s^{-1}')
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Pulse amplitude (a.u.)','FontWeight','bold')
set(gca,'FontSize',18)


%% WT (18 umol) detrended data

s = s_ML;

pkmatun_wt = [];

for i = 2:20
    [pkmatun] = findpeaksmethod_rpoD4(s{i,2},'hrst','MYs_det','lengthMicrons_smooth');
    
    pkmatun_wt = cat(1,pkmatun_wt,pkmatun);
end

%---- Peak amplitude of non-detrended data -----%
Amp = pkmatun_wt(:,5); % Peak amplitude of detrended data
CT = pkmatun_wt(:,4); % Circadian time

clear xs1 yperx1 xs2 meanyperx1 stdyperx1 meanyperx2 stdyperx2 s2

xs1 = unique(CT);
yperx1 = cell(length(xs1),1);
lperx1 = cell(length(xs1),1);

 for j = 1:length(CT)
        whereinxs = find(xs1 == CT(j));
        yperx1{whereinxs} = cat(1,yperx1{whereinxs},Amp(j));
 end

for i = 1:length(yperx1)
    meanyperx1(i) = mean(yperx1{i});
    stdyperx1(i) = std(yperx1{i});
    serrmean1(i) = stdyperx1(i)./sqrt(length(yperx1{i}));
end

xs2 = cat(1,xs1,xs1+24);
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
plot(xs2, meanyperx2, '-','color','c','LineWidth', 2);
set(gca,'XTick',12*(0:40));
grid on
nightmodeon(s_ML{2,2},12,1);
xlim([0 48])
ylim([100 500])
title('18 \mumol m^{-2} s^{-1}')
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Pulse amplitude (a.u.)','FontWeight','bold')
set(gca,'FontSize',18)

%% WT (18 umol) detrended and normalised data
s = s_ML;

pkmatun_wt = [];

for i = 2:20
    [pkmatun] = findpeaksmethod_rpoD4(s{i,2},'hrst','MYs_det','lengthMicrons_smooth');
    
    pkmatun_wt = cat(1,pkmatun_wt,pkmatun);
end


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

meany_ML = mean(meanyperx1);

xs2 = cat(1,xs1,xs1+24);
meanyperx2 = [meanyperx1,meanyperx1];
stdyperx2 = [stdyperx1,stdyperx1];

figure;
curve1 = (meanyperx2 + stdyperx2)./meany_ML;
curve2 = (meanyperx2 - stdyperx2)./meany_ML;
x2 = [xs2', fliplr(xs2')];
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween,'c','LineStyle','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',0.2)
hold on;
plot(xs2, meanyperx2./meany_ML, '-','color','c','LineWidth', 2);
set(gca,'XTick',12*(0:40));
grid on
nightmodeon(s_ML{2,2},12,1);
xlim([0 48])
ylim([0 2.5])
title('18 \mumol m^{-2} s^{-1}')
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Norm. pulse amplitude','FontWeight','bold')
set(gca,'FontSize',18)

%% WT (25 umol) non-detrended data

s = s_HL;

pkmatun_wt = [];

for i = 2:9
    [pkmatun] = findpeaksmethod_rpoD4(s{i,2},'hrst','MYs','lengthMicrons_smooth');
    
    pkmatun_wt = cat(1,pkmatun_wt,pkmatun);
end

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

xs2 = cat(1,xs1,xs1+24);
meanyperx2 = [meanyperx1,meanyperx1];
stdyperx2 = [stdyperx1,stdyperx1];

figure;
curve1 = meanyperx2 + stdyperx2;
curve2 = meanyperx2 - stdyperx2;
x2 = [xs2', fliplr(xs2')];
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween,'c','LineStyle','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',0.2)
hold on;
plot(xs2, meanyperx2, '-','color','c','LineWidth', 2);
set(gca,'XTick',12*(0:40));
grid on
nightmodeon(s_HL{2,2},12,1);
xlim([0 48])
ylim([600 1100])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Pulse amplitude (a.u.)','FontWeight','bold')
title('25 \mumol m^{-2} s^{-1}')
set(gca,'FontSize',18)

%% WT (25 umol) detrended data

s = s_HL;

pkmatun_wt = [];

for i = 2:9
    [pkmatun] = findpeaksmethod_rpoD4(s{i,2},'hrst','MYs_det','lengthMicrons_smooth');
    
    pkmatun_wt = cat(1,pkmatun_wt,pkmatun);
end

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

xs2 = cat(1,xs1,xs1+24);
meanyperx2 = [meanyperx1,meanyperx1];
stdyperx2 = [stdyperx1,stdyperx1];

figure;
curve1 = meanyperx2 + stdyperx2;
curve2 = meanyperx2 - stdyperx2;
x2 = [xs2', fliplr(xs2')];
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween,'c','LineStyle','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',0.2)
hold on;
plot(xs2, meanyperx2, '-','color','c','LineWidth', 2);
set(gca,'XTick',12*(0:40));
grid on
nightmodeon(s_HL{2,2},12,1);
xlim([0 48])
ylim([100 500])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Pulse amplitude (a.u.)','FontWeight','bold')
title('25 \mumol m^{-2} s^{-1}')
set(gca,'FontSize',18)

%% %% WT (25 umol) detrended and normalised data
s = s_HL;

pkmatun_wt = [];

for i = 2:9
    [pkmatun] = findpeaksmethod_rpoD4(s{i,2},'hrst','MYs_det','lengthMicrons_smooth');
    
    pkmatun_wt = cat(1,pkmatun_wt,pkmatun);
end

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

meany_HL = mean(meanyperx1);

xs2 = cat(1,xs1,xs1+24);
meanyperx2 = [meanyperx1,meanyperx1];
stdyperx2 = [stdyperx1,stdyperx1];

figure;
curve1 = (meanyperx2 + stdyperx2)./meany_HL;
curve2 = (meanyperx2 - stdyperx2)./meany_HL;
x2 = [xs2', fliplr(xs2')];
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween,'c','LineStyle','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',0.2)
hold on;
plot(xs2, meanyperx2./meany_HL, '-','color','c','LineWidth', 2);
set(gca,'XTick',12*(0:40));
grid on
nightmodeon(s_HL{2,2},12,1);
xlim([0 48])
ylim([0 2.5])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Norm. pulse amplitude','FontWeight','bold')
title('25 \mumol m^{-2} s^{-1}')
set(gca,'FontSize',18)


%% Fig.S2D WT
s = s_ML;

pkmatun_wt = [];

for i = 2:20
    [pkmatun] = findpeaksmethod_rpoD4(s{i,2},'hrst','MYs_det','lengthMicrons_smooth');
    
    pkmatun_wt = cat(1,pkmatun_wt,pkmatun);
end

%---- Peak width of detrended data -----%
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
h = fill(x2, inBetween,'c','LineStyle','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',0.2)
hold on;
plot(xs2, meanyperx2, '-','color','c','LineWidth', 2);
set(gca,'XTick',12*(0:40));
grid on
nightmodeon(s_ML{2,2},12,1);
xlim([0 48])
ylim([1 7])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Pulse width (h)','FontWeight','bold')
title('WT')
set(gca,'FontSize',18)

%% Fig.S2D delkaiBC
s = s_delkaiBC;

pkmatun_delkaiBC = [];

for i = 2:19
    [pkmatun] = findpeaksmethod_rpoD4(s{i,2},'hrst','MYs_det','lengthMicrons_smooth');
    
    pkmatun_delkaiBC = cat(1,pkmatun_delkaiBC,pkmatun);
end
%---- Peak width of detrended data -----%
Width = pkmatun_delkaiBC(:,6); % Peak width of detrended data
CT = pkmatun_delkaiBC(:,4); % Circadian time

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
h = fill(x2, inBetween,'r','LineStyle','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',0.2)
hold on;
plot(xs2, meanyperx2, '-','color','r','LineWidth', 2);
set(gca,'XTick',12*(0:40));
grid on
nightmodeon(s_delkaiBC{2,2},12,1);
xlim([0 48])
ylim([1 7])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Pulse width (h)','FontWeight','bold')
title('\it\DeltakaiBC')
set(gca,'FontSize',18)

%% Fig.S2E

s = s_kaiBCR393C;

sall_R393C = [];

for k = 2:size(s,1)
s1 = s{k,2};
sall_R393C = cat(2,sall_R393C,s1);
end



Y_wt = [];
L_wt = [];
T_wt = [];
salt_wt = [];

for i = 2:size(s,1)
    [Y salt] = schnitzTables(s{i,2},'MYs_smooth');
    [L salt] = schnitzTables(s{i,2},'lengthMicrons_smooth');
    [T salt] = schnitzTables(s{i,2},'hrst');
    Y_wt = catpad(2,Y_wt,Y);
    L_wt = catpad(2,L_wt,L);
    T_wt = catpad(2,T_wt,T);
    salt_wt = catpad(2,salt_wt,salt);
end


cutoff = 0;
tcutoff = 150;

Y_det = [];

for k = 1:size(T_wt,2)
    t = T_wt(:,k);
    y = Y_wt(:,k);
    idx = intersect(find(t>=cutoff),find(t<=tcutoff));
    y = y(idx);
    t = t(idx);
    f = polyfit(t,y,2);
    ypol = y-t.^2*f(1)-t*f(2)-f(3);
    
    Y_det = catpad(2,Y_det,ypol);
end


figure; 

for i = 1:size(Y_det,2)
    p = patchline(T_wt(:,i),Y_det(:,i),'linestyle','-','edgecolor','[0.1 0.1 0.1]','linewidth',0.5,'edgealpha',0.1);
    p.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold on
end

hold on;


T = [];
Y = [];

for i = 1:size(T_wt,2)
    T = cat(1,T,rmmissing(T_wt(:,i)));
    Y = cat(1,Y,rmmissing(Y_det(:,i)));
end

clear xs1 yperx1 xs2 meanyperx1 stdyperx1 meanyperx2 stdyperx2 s2

xs1 = unique(T);
yperx1 = cell(length(xs1),1);
lperx1 = cell(length(xs1),1);

 for j = 1:length(T)
        whereinxs = find(xs1 == T(j));
        yperx1{whereinxs} = cat(1,yperx1{whereinxs},Y(j));
 end

 
for i = 1:length(yperx1)
    meanyperx1(i) = mean(yperx1{i});
    stdyperx1(i) = std(yperx1{i});
    serrmean1(i) = stdyperx1(i)./sqrt(length(yperx1{i}));
end


curve1 = meanyperx1 + stdyperx1;
curve2 = meanyperx1 - stdyperx1;
x2 = [xs1', fliplr(xs1')];
inBetween = [curve1, fliplr(curve2)];

h = fill(x2, inBetween,'r','LineStyle','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',0.2)
hold on;

plot(xs1, meanyperx1, '-','color','r','LineWidth', 2);

set(gca,'XTick',7.5*(0:40));
grid on
nightmodeon15(s_ML{20,2},7.5,1);
xlim([0 45])
ylim([-400 600])
xlabel('Time in LL (h)','FontWeight','bold')
ylabel('Dtr. P_{\itkaiBC}-EYFP (a.u.)','FontWeight','bold')
title('KaiCR393C')
set(gca,'FontSize',14)