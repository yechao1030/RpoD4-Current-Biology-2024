%% Fig.3A
clear
load('Fig3A_data_1.mat')

s = s_VLL;
sall_VLL = [];

for k = 2:11
s1 = s{k,2};
sall_VLL = cat(2,sall_VLL,s1);
end



Y_VLL = [];
L_VLL = [];
T_VLL = [];
salt_VLL = [];

for i = 2:11
    [Y salt] = schnitzTables(s{i,2},'MYs_det');
    [L salt] = schnitzTables(s{i,2},'lengthMicrons_smooth');
    [T salt] = schnitzTables(s{i,2},'hrst');
    Y_VLL = catpad(2,Y_VLL,Y);
    L_VLL = catpad(2,L_VLL,L);
    T_VLL = catpad(2,T_VLL,T);
    salt_VLL = catpad(2,salt_VLL,salt);
end


s = s_ML;
sall_ML = [];

for k = 2:size(s,1)
s1 = s{k,2};
sall_ML = cat(2,sall_ML,s1);
end



Y_ML = [];
L_ML = [];
T_ML = [];
salt_ML = [];

for i = 2:size(s,1)
    [Y salt] = schnitzTables(s{i,2},'MYs_det');
    [L salt] = schnitzTables(s{i,2},'lengthMicrons_smooth');
    [T salt] = schnitzTables(s{i,2},'hrst');
    Y_ML = catpad(2,Y_ML,Y);
    L_ML = catpad(2,L_ML,L);
    T_ML = catpad(2,T_ML,T);
    salt_ML = catpad(2,salt_ML,salt);
end


s = s_HL;
sall_HL = [];

for k = 2:size(s,1)
s1 = s{k,2};
sall_HL = cat(2,sall_HL,s1);
end



Y_HL = [];
L_HL = [];
T_HL = [];
salt_HL = [];

for i = 2:size(s,1)
    [Y salt] = schnitzTables(s{i,2},'MYs_det');
    [L salt] = schnitzTables(s{i,2},'lengthMicrons_smooth');
    [T salt] = schnitzTables(s{i,2},'hrst');
    Y_HL = catpad(2,Y_HL,Y);
    L_HL = catpad(2,L_HL,L);
    T_HL = catpad(2,T_HL,T);
    salt_HL = catpad(2,salt_HL,salt);
end

%%% Plot Figure 3A panels %%%

figure;

y = rmmissing(Y_VLL(:,100));
x = rmmissing(T_VLL(:,100));
[pk,lk] = findpeaks(y,x,'MinPeakProminence',25);
plot(x,y/meany_VLL,'-','Color',[0.00,0.45,0.74],'LineWidth',2);
hold on
plot(lk,pk./meany_VLL,'o','Color',[0.00,0.45,0.74],'LineWidth',2, 'MarkerSize',8);
xlim([0 72])
ylim([-0.5 2])
set(gca,'XTick',12*(0:40));
grid on
title('10 \muE m^{-2} s^{-1}')
hold off


figure;

y = rmmissing(Y_ML(:,300));
x = rmmissing(T_ML(:,300));
[pk,lk] = findpeaks(y,x,'MinPeakProminence',100);
plot(x,y/meany_ML,'-','Color',[0.47,0.67,0.19],'LineWidth',2);
hold on
plot(lk,pk./meany_ML,'o','Color',[0.47,0.67,0.19],'LineWidth',2, 'MarkerSize',8);

xlim([0 72])
ylim([-0.5 2])
set(gca,'XTick',12*(0:40));
grid on
ylabel('Norm. detr. P_{\itrpoD4}-EYFP','FontWeight','bold')
title('18 \muE m^{-2} s^{-1}')


figure;

y = rmmissing(Y_HL(:,409));
x = rmmissing(T_HL(:,409));
[pk,lk] = findpeaks(y,x,'MinPeakProminence',100);
plot(x,y/meany_HL,'-','Color',[0.85,0.33,0.10],'LineWidth',2);
hold on
plot(lk,pk./meany_HL,'o','Color',[0.85,0.33,0.10],'LineWidth',2, 'MarkerSize',8);
xlim([0 72])
ylim([-0.5 2])
set(gca,'XTick',12*(0:40));
grid on
xlabel('Time in LL (h)','fontweight','bold')
title('25 \muE m^{-2} s^{-1}')


%% Fig.3B 
%%% Very low light %%%

load('Fig3A_data_1.mat')

s = s_VLL;
sall_wt = [];

for k = 2:11
s1 = s{k,2};
sall_wt = cat(2,sall_wt,s1);
end



Y_wt = [];
L_wt = [];
T_wt = [];
salt_wt = [];

for i = 2:11
    [Y salt] = schnitzTables(s{i,2},'MYs_det');
    [L salt] = schnitzTables(s{i,2},'lengthMicrons_smooth');
    [T salt] = schnitzTables(s{i,2},'hrst');
    Y_wt = catpad(2,Y_wt,Y);
    L_wt = catpad(2,L_wt,L);
    T_wt = catpad(2,T_wt,T);
    salt_wt = catpad(2,salt_wt,salt);
end


figure
CT = (0:0.75:23.75);
%cmp = parula(34);
cmp = parula(16);
cmp = cat(1,flipud(cmp),cmp);
for k = 1:size(T_wt,2)
    y = rmmissing(Y_wt(:,k));
    t = rmmissing(T_wt(:,k));
    [pks locs w p] = findpeaks(y,t,'MinPeakProminence',25);
    [pks0 idx0] = findpeaks(y,'MinPeakProminence',25);
    for n = 1:length(locs)
        ct = mod(locs(n),24);
        idx = find(CT==ct);
        plot(t([max(idx0(n)-9,1):min(idx0(n)+9,length(y))]),y([max(idx0(n)-9,1):min(idx0(n)+9,length(y))])./meany_VLL,'color',[cmp(idx,:)])% normalise by mean pulse height
        
        
        hold on
    end
end
set(gca,'XTick',12*(0:40));
grid on;
xlim([0 120])
ylim([-0.5 2.5])
colormap(cmp)

colorbar('Ticks',[0:1/32:1],...
    'TickLabels',{'0','','','','3','','','','6','','','','9','','',...
    '','12','','','','15','','','','18','','','','21','','','','24'})

xlabel('Time in LL (h)','FontWeight','bold')
ylabel('Norm. detr. P_{\itrpoD4}-EYFP pulse','FontWeight','bold')
title('10 \mumol m^{-2} s^{-1} (ccd = 22.5 h)')



%%% Medium light WT %%%
s = s_ML;
sall_wt = [];

for k = 2:size(s,1)
s1 = s{k,2};
sall_wt = cat(2,sall_wt,s1);
end



Y_wt = [];
L_wt = [];
T_wt = [];
salt_wt = [];

for i = 2:size(s,1)
    [Y salt] = schnitzTables(s{i,2},'MYs_det');
    [L salt] = schnitzTables(s{i,2},'lengthMicrons_smooth');
    [T salt] = schnitzTables(s{i,2},'hrst');
    Y_wt = catpad(2,Y_wt,Y);
    L_wt = catpad(2,L_wt,L);
    T_wt = catpad(2,T_wt,T);
    salt_wt = catpad(2,salt_wt,salt);
end

figure
CT = (0:0.75:23.75);
cmp = parula(16);
cmp = cat(1,flipud(cmp),cmp);
for k = 1:size(T_wt,2)
    y = rmmissing(Y_wt(:,k));
    t = rmmissing(T_wt(:,k));
    [pks locs w p] = findpeaks(y,t,'MinPeakProminence',100);
    [pks0 idx0] = findpeaks(y,'MinPeakProminence',100);


    for n = 1:length(locs)
        ct = mod(locs(n),24);
        idx = find(CT==ct);
        plot(t([max(idx0(n)-9,1):min(idx0(n)+9,length(y))]),y([max(idx0(n)-9,1):min(idx0(n)+9,length(y))])./meany_ML,'color',[cmp(idx,:)])% normalise by mean pulse height
        
        
        hold on
    end
end
set(gca,'XTick',12*(0:40));
grid on;
xlim([0 96])
ylim([-1 2.5])
colormap(cmp)
colorbar('Ticks',[0:1/32:1],...
    'TickLabels',{'0','','','','3','','','','6','','','','9','','',...
    '','12','','','','15','','','','18','','','','21','','','','24'})
xlabel('Time in LL (h)','FontWeight','bold')
ylabel('Norm. detr. P_{\itrpoD4}-EYFP pulse','FontWeight','bold')
title('18 \mumol m^{-2} s^{-1} (ccd = 14.25 h)')



%%% High light %%%
s = s_HL;
sall_wt = [];

for k = 2:size(s,1)
s1 = s{k,2};
sall_wt = cat(2,sall_wt,s1);
end



Y_wt = [];
L_wt = [];
T_wt = [];
salt_wt = [];

for i = 2:size(s,1)
    [Y salt] = schnitzTables(s{i,2},'MYs_det');
    [L salt] = schnitzTables(s{i,2},'lengthMicrons_smooth');
    [T salt] = schnitzTables(s{i,2},'hrst');
    Y_wt = catpad(2,Y_wt,Y);
    L_wt = catpad(2,L_wt,L);
    T_wt = catpad(2,T_wt,T);
    salt_wt = catpad(2,salt_wt,salt);
end

figure
CT = (0:0.75:23.75);
cmp = parula(16);
cmp = cat(1,flipud(cmp),cmp);
for k = 1:size(T_wt,2)
    y = rmmissing(Y_wt(:,k));
    t = rmmissing(T_wt(:,k));
    [pks locs w p] = findpeaks(y,t,'MinPeakProminence',100);
    [pks0 idx0] = findpeaks(y,'MinPeakProminence',100);
    for n = 1:length(locs)
        ct = mod(locs(n),24);
        idx = find(CT==ct);
        plot(t([max(idx0(n)-9,1):min(idx0(n)+9,length(y))]),y([max(idx0(n)-9,1):min(idx0(n)+9,length(y))])./meany_HL,'color',[cmp(idx,:)])%normalised by mean pulse height
        
        hold on
    end
end
set(gca,'XTick',12*(0:40));
grid on;
xlim([0 72])
ylim([-1 2.5])
colormap(cmp)
colorbar('Ticks',[0:1/32:1],...
    'TickLabels',{'0','','','','3','','','','6','','','','9','','',...
    '','12','','','','15','','','','18','','','','21','','','','24'})
xlabel('Time in LL (h)','FontWeight','bold')
ylabel('Norm. detr. P_{\itrpoD4}-EYFP pulse','FontWeight','bold')
title('25 \mumol m^{-2} s^{-1} (ccd = 10.5 h)')

%% Fig.3C
clear
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


%% Fig.3D

clear 

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

%%% Find rpoD4 expression peaks in the KaiCR393C mutant %%

pkmatun_wt = [];

for i = 2:size(s,1)
    [pkmatun] = findpeaksmethod_rpoD4(s{i,2},'hrst','MYs_det','lengthMicrons_smooth',15,0,30,100);
    
    pkmatun_wt = cat(1,pkmatun_wt,pkmatun);
end


%%% Plot Fig3D (no binning normalised) %%%
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
