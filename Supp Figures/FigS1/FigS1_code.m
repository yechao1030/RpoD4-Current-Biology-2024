%% Fig.S1A (delkaiBC)
clear
load('schnitzes-2021-07-20-V1_corrected.mat')
[Y_delkaiBC salt_delkaiBC] = schnitzTables(schnitz1{18,2},'MYs');
[L_delkaiBC salt_delkaiBC] = schnitzTables(schnitz1{18,2},'lengthMicrons_smooth');
[T_delkaiBC salt_delkaiBC] = schnitzTables(schnitz1{18,2},'hrst');

figure;

subplot(2,1,1);

plot(T_delkaiBC,Y_delkaiBC,'k-');
hold on;
[pk,lk] = findpeaks((Y_delkaiBC(:,15)),T_delkaiBC(:,15),'MinPeakProminence',100);
plot(T_delkaiBC(:,15),Y_delkaiBC(:,15),'r-',lk,pk,'yo','LineWidth',2, 'MarkerSize',8);

nightmodeon(schnitz1{18,2},12,1);
xlabel('')
ylabel('P_{rpoD4}-EYFP (a.u.)','fontweight','bold')
title('\it\DeltakaiBC')
ylim([300 1100])
xlim([24 96])
set(gca,'XTick',12*(0:20));
set(gca,'YTick',100*(0:20));
set(gca,'FontSize',12);
grid on;

subplot(2,1,2);
plot(T_delkaiBC,L_delkaiBC,'k-');
hold on
TF = islocalmin(L_delkaiBC(:,15));
x = T_delkaiBC(:,15);y = L_delkaiBC(:,15);
plot(x,y,'r-',x(TF),y(TF),'yo','LineWidth',2, 'MarkerSize',8);
nightmodeon(schnitz1{18,2},12,1);
xlabel('Time in LL (h)','fontweight','bold')
ylabel('Cell length (\mum)','fontweight','bold')
ylim([2 10])
xlim([24 96])
set(gca,'XTick',12*(0:20));
set(gca,'FontSize',12);
grid on;

%% Fig.S1C (WT)
clear

load('FigS1_data.mat')

[Y salt] = schnitzTables(s_wt{8,2},'MYs');
[T salt] = schnitzTables(s_wt{8,2},'hrst');
[L salt] = schnitzTables(s_wt{8,2},'lengthMicrons');

figure;

subplot(2,1,1);
plot(T,Y,'k-');
hold on;
[pk,lk] = findpeaks((Y(:,32)),T(:,32),'MinPeakProminence',40);
plot(T(:,32),Y(:,32),'c-',lk,pk,'mo','LineWidth',2, 'MarkerSize',8);
nightmodeon(s_wt{16,2},12,1);
xlabel('')
ylabel('RpoD4-EYFP (a.u.)','fontweight','bold')
title('WT')
ylim([200 900])
xlim([0 72])
set(gca,'XTick',12*(0:20));
set(gca,'YTick',100*(0:20));
set(gca,'FontSize',12);
grid on;

subplot(2,1,2);
plot(T,L,'k-');
hold on
TF = islocalmin(L(:,32));
x = T(:,32);y = L(:,32);
plot(x,y,'c-',x(TF),y(TF),'mo','LineWidth',2, 'MarkerSize',8);
nightmodeon(s_wt{16,2},12,1);
xlabel('Time in LL (h)','fontweight','bold')
ylabel('Cell length (\mum)','fontweight','bold')
ylim([2 10])
xlim([0 72])
set(gca,'XTick',12*(0:20));
set(gca,'FontSize',12);
grid on;


%% Fig.S1D (delkaiBC)
[Y salt] = schnitzTables(s_delkaiBC{13,2},'MYs');
[L salt] = schnitzTables(s_delkaiBC{13,2},'lengthMicrons_smooth');
[T salt] = schnitzTables(s_delkaiBC{13,2},'hrst');

figure;

subplot(2,1,1);

plot(T,Y,'k-');
hold on;
[pk,lk] = findpeaks((Y(:,119)),T(:,119),'MinPeakProminence',40);
plot(T(:,119),Y(:,119),'r-',lk,pk,'yo','LineWidth',2, 'MarkerSize',8);
nightmodeon(s_delkaiBC{13,2},12,1);
xlabel('')
ylabel('RpoD4-EYFP (a.u.)','fontweight','bold')
title('\it\DeltakaiBC')
ylim([200 900])
xlim([0 72])
set(gca,'XTick',12*(0:20));
set(gca,'YTick',100*(0:20));
set(gca,'FontSize',12);
grid on;

subplot(2,1,2);
plot(T,L,'k-');
hold on
TF = islocalmin(L(:,119));
x = T(:,119);y = L(:,119);
plot(x,y,'r-',x(TF),y(TF),'yo','LineWidth',2, 'MarkerSize',8);
nightmodeon(s_delkaiBC{13,2},12,1);
xlabel('Time in LL (h)','fontweight','bold')
ylabel('Cell length (\mum)','fontweight','bold')
ylim([2 10])
xlim([0 72])
set(gca,'XTick',12*(0:20));
set(gca,'FontSize',12);
grid on;


%% Fig.S1E (WT PrpoD4-EYFP vs ToD)
clear

load('FigS1_data.mat')

T_wt = [];
T_wt_all = [];
Y_wt = [];
Y_wt_all = [];

for j = 2:length(s_wt)
    for i = 1:length(s_wt{j,2})
        if s_wt{j,2}(i).E ~= 0 && s_wt{j,2}(i).P ~= 0
            A = s_wt{j,2}(i).E;
            B = s_wt{j,2}(i).D;
            t1 = [s_wt{j,2}(i).hrst-max(s_wt{j,2}(i).hrst)-0.75,s_wt{j,2}(A).hrst-min(s_wt{j,2}(A).hrst)];
            y1 = [s_wt{j,2}(i).MYs_det,s_wt{j,2}(A).MYs_det];
            
            
            t2 = [s_wt{j,2}(i).hrst-max(s_wt{j,2}(i).hrst)-0.75,s_wt{j,2}(B).hrst-min(s_wt{j,2}(B).hrst)];
            y2 = [s_wt{j,2}(i).MYs_det,s_wt{j,2}(B).MYs_det];
                
           
            T_wt = catpad(2,T_wt,t1');
            T_wt = catpad(2,T_wt,t2');
            Y_wt = catpad(2,Y_wt,y1');
            Y_wt = catpad(2,Y_wt,y2');

            T_wt_all = cat(1,T_wt_all,t1');
            T_wt_all = cat(1,T_wt_all,t2');
            Y_wt_all = cat(1,Y_wt_all,y1');
            Y_wt_all = cat(1,Y_wt_all,y2');
        end
    end
end

xs = unique(T_wt_all);

yperx = cell(length(xs),1);
meanyperx = nan(1,length(yperx));
stdyperx = nan(1,length(yperx));
serrmeany = nan(1,length(yperx));


for j = 1:length(T_wt_all)
        whereinxs = find(xs == T_wt_all(j));
        yperx{whereinxs} = cat(1,yperx{whereinxs},Y_wt_all(j));
end


for i = 1:length(yperx)
    meanyperx(i) = mean(yperx{i});
    stdyperx(i) = std(yperx{i});
    serrmeany(i) = stdyperx(i)./sqrt(length(yperx{i}));
end



figure;
plot(T_wt,Y_wt,'.-');
hold on;

curve1 = meanyperx + stdyperx;
curve2 = meanyperx - stdyperx;
x2 = [xs', fliplr(xs')];
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween,'y','LineStyle','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',0.2)
hold on;

plot(xs,meanyperx,'y-','LineWidth',4)


title('WT')
xlabel('Time after cell division (h)','FontWeight','bold');
ylabel('Detrended RpoD4-EYFP (a.u.)','FontWeight','bold')
grid on;
xline(0,'k--','LineWidth',4)
xlim([-9,9]);
ylim([-200,400]);
set(gca,'XTick',3*(-40:40));
set(gca,'FontSize',14)

%% Fig.1F (delkaiBC PrpoD4-EYFP vs ToD)
clear

load('FigS1_data.mat')

T_delkaiBC = [];
T_delkaiBC_all = [];
Y_delkaiBC = [];
Y_delkaiBC_all = [];

for j = 2:length(s_delkaiBC)
    for i = 1:length(s_delkaiBC{j,2})
        if s_delkaiBC{j,2}(i).E ~= 0 && s_delkaiBC{j,2}(i).P ~= 0
            A = s_delkaiBC{j,2}(i).E;
            B = s_delkaiBC{j,2}(i).D;
            t1 = [s_delkaiBC{j,2}(i).hrst-max(s_delkaiBC{j,2}(i).hrst)-0.75,s_delkaiBC{j,2}(A).hrst-min(s_delkaiBC{j,2}(A).hrst)];
            y1 = [s_delkaiBC{j,2}(i).MYs_det,s_delkaiBC{j,2}(A).MYs_det];
            
            t2 = [s_delkaiBC{j,2}(i).hrst-max(s_delkaiBC{j,2}(i).hrst)-0.75,s_delkaiBC{j,2}(B).hrst-min(s_delkaiBC{j,2}(B).hrst)];
            y2 = [s_delkaiBC{j,2}(i).MYs_det,s_delkaiBC{j,2}(B).MYs_det];          
           
            T_delkaiBC = catpad(2,T_delkaiBC,t1');
            T_delkaiBC = catpad(2,T_delkaiBC,t2');
            Y_delkaiBC = catpad(2,Y_delkaiBC,y1');
            Y_delkaiBC = catpad(2,Y_delkaiBC,y2');

            T_delkaiBC_all = cat(1,T_delkaiBC_all,t1');
            T_delkaiBC_all = cat(1,T_delkaiBC_all,t2');
            Y_delkaiBC_all = cat(1,Y_delkaiBC_all,y1');
            Y_delkaiBC_all = cat(1,Y_delkaiBC_all,y2');
        end
    end
end

xs = unique(T_delkaiBC_all);

yperx = cell(length(xs),1);
meanyperx = nan(1,length(yperx));
stdyperx = nan(1,length(yperx));
serrmeany = nan(1,length(yperx));


for j = 1:length(T_delkaiBC_all)
        whereinxs = find(xs == T_delkaiBC_all(j));
        yperx{whereinxs} = cat(1,yperx{whereinxs},Y_delkaiBC_all(j));
end


for i = 1:length(yperx)
    meanyperx(i) = mean(yperx{i});
    stdyperx(i) = std(yperx{i});
    serrmeany(i) = stdyperx(i)./sqrt(length(yperx{i}));
end



figure;
plot(T_delkaiBC,Y_delkaiBC,'.-');
hold on;

curve1 = meanyperx + stdyperx;
curve2 = meanyperx - stdyperx;
x2 = [xs', fliplr(xs')];
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween,'y','LineStyle','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',0.2)
hold on;

plot(xs,meanyperx,'y-','LineWidth',4)

title('\it\DeltakaiBC')
xlabel('Time after cell division (h)','FontWeight','bold');
ylabel('Detrended RpoD4-EYFP (a.u.)','FontWeight','bold')
grid on;
xline(0,'k--','LineWidth',4)
xlim([-9,9]);
ylim([-200,400]);
set(gca,'XTick',3*(-40:40));
set(gca,'FontSize',14)


%% WT Detrend rpoD4 signal using basal expression level for each movie
clear

load('FigS1_data.mat')

for i = 2:16
    xmins = []; 
    tmins = [];
    for j = 1:length(s_wt{i,2})
        x = s_wt{i,2}(j).MYs; 
        minx = min(x); 
        indx = find(x == minx); 
        xmins = cat(2,xmins,minx);
        tmins = cat(2,tmins,s_wt{i,2}(j).hrst(indx));
    end

    tmax = max([s_wt{i,2}.hrst]);

    xmins(tmins > tmax - 6) = [];
    tmins(tmins > tmax - 6) = [];

f = polyfit(tmins,xmins,2);


    for j = 1:length(s_wt{i,2})
        ts = s_wt{i,2}(j).hrst; 
        ys = s_wt{i,2}(j).MYs; 
        s_wt{i,2}(j).MYs_det = ys - ts.^2*f(1) - ts*f(2) - f(3);
    end
end

%% delkaiBC Detrend rpoD4 signal using basal expression level for each movie
for i = 2:17
    xmins = []; 
    tmins = [];
    for j = 1:length(s_delkaiBC{i,2})
        x = s_delkaiBC{i,2}(j).MYs; 
        minx = min(x); 
        indx = find(x == minx); 
        xmins = cat(2,xmins,minx);
        tmins = cat(2,tmins,s_delkaiBC{i,2}(j).hrst(indx));
    end

    tmax = max([s_delkaiBC{i,2}.hrst]);

    xmins(tmins > tmax - 6) = [];
    tmins(tmins > tmax - 6) = [];

f = polyfit(tmins,xmins,2);


    for j = 1:length(s_delkaiBC{i,2})
        ts = s_delkaiBC{i,2}(j).hrst; 
        ys = s_delkaiBC{i,2}(j).MYs; 
        s_delkaiBC{i,2}(j).MYs_det = ys - ts.^2*f(1) - ts*f(2) - f(3);
    end
end
%%
pkmatun_wt = [];


for i = 2:length(s_wt)
    [pkmatun] = findpeaksmethod_RpoD4_v1(s_wt{i,2},'hrst','MYs_det','lengthMicrons_smooth');
    
    pkmatun_wt = cat(1,pkmatun_wt,pkmatun);
end

%%
pkmatun_delkaiBC = [];

for i = 2:length(s_delkaiBC)
    [pkmatun] = findpeaksmethod_RpoD4_v1(s_delkaiBC{i,2},'hrst','MYs_det','lengthMicrons_smooth');
    
    pkmatun_delkaiBC = cat(1,pkmatun_delkaiBC,pkmatun);
end


%% Fig.S1G (WT)
%---- Peak height of detrended data -----%
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
    ncells(i) = length(yperx1{i});
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
nightmodeon(s_wt{16,2},12,1);
xlim([0 48])
ylim([0 250])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Pulse amplitude (a.u.)','FontWeight','bold')
title('WT')
set(gca,'FontSize',18)


%% Fig.S1H (WT)
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
nightmodeon(s_wt{16,2},12,1);
xlim([0 48])
ylim([0 6])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Pulse width (h)','FontWeight','bold')
set(gca,'FontSize',18)

%% Fig.S1I (delkaiBC)
%---- Peak height of detrended data -----%
Amp = pkmatun_delkaiBC(:,5);
CT = pkmatun_delkaiBC(:,4);

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
    ncells(i) = length(yperx1{i});
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
nightmodeon(s_delkaiBC{13,2},12,1);
xlim([0 48])
ylim([0 250])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Pulse amplitude (a.u.)','FontWeight','bold')
title('\it\DeltakaiBC')
set(gca,'FontSize',18)

%% Fig.S1J (delkaiBC)
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
nightmodeon(s_delkaiBC{13,2},12,1);
xlim([0 48])
ylim([0 6])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Pulse width (h)','FontWeight','bold')
title('\it\DeltakaiBC')
set(gca,'FontSize',18)
