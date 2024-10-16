%% Fig.4B
clear 
load("Fig4AB_data.mat")

%%% WT %%%
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

%%% delrpoD4 %%%
for i = 2:19
    xmins = []; 
    tmins = [];
    for j = 1:length(s_delrpoD4{i,2})
        x = s_delrpoD4{i,2}(j).MYs; 
        minx = min(x); 
        indx = find(x == minx); 
        xmins = cat(2,xmins,minx);
        tmins = cat(2,tmins,s_delrpoD4{i,2}(j).hrst(indx));
    end

    tmax = max([s_delrpoD4{i,2}.hrst]);

    xmins(tmins > tmax - 6) = [];
    tmins(tmins > tmax - 6) = [];

f = polyfit(tmins,xmins,2);


    for j = 1:length(s_delrpoD4{i,2})
        ts = s_delrpoD4{i,2}(j).hrst; 
        ys = s_delrpoD4{i,2}(j).MYs; 
        s_delrpoD4{i,2}(j).MYs_det = ys - ts.^2*f(1) - ts*f(2) - f(3);
    end
end

%%% Get rpoD4 peak information %%%
%%% WT %%%
s = s_ML;

pkmatun_wt = [];

for i = 2:20
    [pkmatun] = findpeaksmethod_rpoD4(s{i,2},'hrst','MYs_det','lengthMicrons_smooth',24,0,144,100);
    
    pkmatun_wt = cat(1,pkmatun_wt,pkmatun);
end

%%% delrpoD4 %%%
s = s_delrpoD4;

pkmatun_delrpoD4 = [];

for i = 2:19
    [pkmatun] = findpeaksmethod_rpoD4(s{i,2},'hrst','MYs_det','lengthMicrons_smooth',24,0,144,100);
    
    pkmatun_delrpoD4 = cat(1,pkmatun_delrpoD4,pkmatun);
end



%%% Plot Fig4B %%%
%---- WT Peak amplitude of detrended data -----%
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


%---- delrpoD4 Peak amplitude of detrended data -----%
clear Amp CT

Amp = pkmatun_delrpoD4(:,5); % Peak amplitude of detrended data
CT = pkmatun_delrpoD4(:,4); % Circadian time

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

xs2 = cat(1,xs1,xs1+24);
meanyperx2 = [meanyperx1,meanyperx1];
stdyperx2 = [stdyperx1,stdyperx1];

curve1 = meanyperx2 + stdyperx2;
curve2 = meanyperx2 - stdyperx2;
x2 = [xs2', fliplr(xs2')];
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween,'b','LineStyle','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',0.2)
hold on;
plot(xs2, meanyperx2, '-','color','b','LineWidth', 2);



set(gca,'XTick',12*(0:40));
grid on
nightmodeon(s_ML{2,2},12,1);
xlim([0 48])
ylim([100 600])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Pulse amplitude (a.u.)','FontWeight','bold')
legend({'WT','\it\DeltarpoD4'},'Location','northeast')
set(gca,'FontSize',18)


%% Fig.4A
%%% Compile data into one Schnitzcell file %%%
sall_wt = [];
sall_delrpoD4 = [];

for k = 2:size(s_ML,1)
s1 = s_ML{k,2};
sall_wt = cat(2,sall_wt,s1);
end

for k = 2:size(s_delrpoD4,1)
s1 = s_delrpoD4{k,2};
sall_delrpoD4 = cat(2,sall_delrpoD4,s1);
end

%%% Get cell properties inclduing cell length at birth (i.e. sb) %%%
[age_wt, ctob_wt,ctod_wt,yob_wt,sb_wt, sd_wt, size_diff_wt, relsize_diff_wt, elongrate_wt, cellno_wt] = schnitzphysprops2(sall_wt,[],'hrst','lengthMicrons_smooth','MYs_det');
[age_delrpoD4, ctob_delrpoD4,ctod_delrpoD4,yob_delrpoD4,sb_delrpoD4, sd_delrpoD4, size_diff_delrpoD4, relsize_diff_delrpoD4, elongrate_delrpoD4, cellno_delrpoD4] = schnitzphysprops2(sall_delrpoD4,[],'hrst','lengthMicrons_smooth','MYs_det');

%%% Plot Fig4A %%%
clear xs1 yperx1 lperx1

xs1 = unique(ctob_wt);
yperx1 = cell(length(xs1),1);
lperx1 = cell(length(xs1),1);

 for j = 1:length(ctob_wt)
        whereinxs = find(xs1 == ctob_wt(j));
        yperx1{whereinxs} = cat(1,yperx1{whereinxs},sb_wt(j));
 end

clear meanyperx1 stdyperx1 serrmean1

for i = 1:length(yperx1)
    meanyperx1(i) = mean(yperx1{i});
    stdyperx1(i) = std(yperx1{i});
    serrmean1(i) = stdyperx1(i)./sqrt(length(yperx1{i}));
end

clear xs2 meanyperx2 stdyperx2

xs2 = cat(2,xs1,xs1+24);
meanyperx2 = [meanyperx1,meanyperx1];
stdyperx2 = [stdyperx1,stdyperx1];

clear curve1 curv2 x2 inBetween h

figure
curve1 = meanyperx2 + stdyperx2;
curve2 = meanyperx2 - stdyperx2;
x2 = [xs2, fliplr(xs2)];
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween,'c','LineStyle','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',0.2)
hold on;
plot(xs2, meanyperx2, '-','color','c','LineWidth', 2);


%---- delrpoD4 Peak amplitude of detrended data -----%

clear xs1 yperx1 lperx1

xs1 = unique(ctob_delrpoD4);
yperx1 = cell(length(xs1),1);
lperx1 = cell(length(xs1),1);

 for j = 1:length(ctob_delrpoD4)
        whereinxs = find(xs1 == ctob_delrpoD4(j));
        yperx1{whereinxs} = cat(1,yperx1{whereinxs},sb_delrpoD4(j));
 end

clear meanyperx1 stdyperx1 serrmean1

for i = 1:length(yperx1)
    meanyperx1(i) = mean(yperx1{i});
    stdyperx1(i) = std(yperx1{i});
    serrmean1(i) = stdyperx1(i)./sqrt(length(yperx1{i}));
end

clear xs2 meanyperx2 stdyperx2

xs2 = cat(2,xs1,xs1+24);
meanyperx2 = [meanyperx1,meanyperx1];
stdyperx2 = [stdyperx1,stdyperx1];

clear curve1 curv2 x2 inBetween h

curve1 = meanyperx2 + stdyperx2;
curve2 = meanyperx2 - stdyperx2;
x2 = [xs2, fliplr(xs2)];
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween,'b','LineStyle','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',0.2)
hold on;
plot(xs2, meanyperx2, '-','color','b','LineWidth', 2);


set(gca,'XTick',12*(0:40));
grid on
nightmodeon(s_ML{2,2},12,1);
xlim([0 48])
ylim([2.5 5])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Birth length (\mum)','FontWeight','bold')
legend({'WT','\it\DeltarpoD4'},'Location','northeast')
set(gca,'FontSize',18)

%% Fig.4C
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


%% Fig.4D
clear
load("Ptrc_sigma_OX_data.mat")

figure;

LB = [A(:,11);A(:,1);A(:,2);A(:,3);A(:,4);A(:,5);A(:,6);A(:,8);A(:,7);A(:,9);...
    A(:,10);A(:,12);A(:,13);A(:,14);A(:,15);A(:,16);A(:,17);A(:,19);A(:,18);...
    A(:,20);A(:,21)];
groupIdx = [ones(size(A(:,11)));2*ones(size(A(:,1)));3*ones(size(A(:,2)));...
    4*ones(size(A(:,3)));5*ones(size(A(:,4)));6*ones(size(A(:,5)));...
    7*ones(size(A(:,6)));8*ones(size(A(:,8)));9*ones(size(A(:,7)));...
    10*ones(size(A(:,9)));11*ones(size(A(:,10)));12*ones(size(A(:,12)));...
    13*ones(size(A(:,13)));14*ones(size(A(:,14)));15*ones(size(A(:,15)));...
    16*ones(size(A(:,16)));17*ones(size(A(:,17)));18*ones(size(A(:,19)));...
    19*ones(size(A(:,18)));20*ones(size(A(:,20)));21*ones(size(A(:,21)))];
violin = violinplot(LB,groupIdx,'ViolinColor',[0.8,0.8,0.8]);

violin(1,14).ViolinColor = [0.9290 0.6940 0.1250];
violin(1,15).ViolinColor = [1,0,0];
violin(1,16).ViolinColor = [0.9290 0.6940 0.1250];
violin(1,18).ViolinColor = [0.9290 0.6940 0.1250];

for i = 1:21
    violin(1,i).ScatterPlot.SizeData = 10;
end


xticklabels({'WT','\itrpoD1','\itrpoD2','\itrpoD3','\itrpoD4','\itrpoD5'...
    ,'\itrpoD6','\itsigF','\itsigF2','\itsigG','\itsigI','\itrpoD1OX'...
    ,'\itrpoD2OX','\itrpoD3OX','\itrpoD4OX','\itrpoD5OX','\itrpoD6OX'...
    ,'\itsigFOX','\itsigF2OX','\itsigGOX','\itsigIOX'});
xtickangle(45)
ylabel('Cell length (\mum)','fontweight','bold')
set(gca,'FontSize',16);

%% Fig.4E
clear
load("Ptrc-rpoD4 dose response curve data.mat")

figure;

[hillCoeff ec50]=doseResponse(d,r);
xlabel('[IPTG] \muM','FontWeight','bold')
ylabel('Cell length (\mum)','FontWeight','bold')
grid on
set(gca,'FontSize',16)

xlim([0.1 10000])