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


%     plotschnitzme(s{i,2},'hrst','MYs',[],'b-'),hold on, plot(tmins,xmins,'*'),hold on
%     ts = 0:80; plot(ts,ts.^2*f(1) + ts*f(2) + f(3),'k','LineWidth',2),hold off
%     pause
    
    for j = 1:length(s{i,2})
        ts = s{i,2}(j).hrst; 
        ys = s{i,2}(j).MYs; 
        s{i,2}(j).MYs_det = ys - ts.^2*f(1) - ts*f(2) - f(3);
    end
end
%%
%s = s;

pkmatun_wt = [];
% pkmatun_delD4 = [];
% pkmatun_delkaiBC = [];

for i = 2:size(s,1)
    [pkmatun] = findpeaksmethod_rpoD4_R393C(s{i,2},'hrst','MYs_det','lengthMicrons_smooth');
    
    pkmatun_wt = cat(1,pkmatun_wt,pkmatun);
end

% for i = 2:size(s_delrpoD4,1)
%     [pkmatun] = findpeaksmethod_rpoD4(s_wt{i,2},'hrst','MYs','MRs','lengthMicrons_smooth');
%     
%     pkmatun_delD4 = cat(1,pkmatun_delD4,pkmatun);
% end
% 
% for i = 2:length(s_delkaiBC)
%     [pkmatun] = findpeaksmethod_rpoD4(s_delkaiBC{i,2},'hrst','MYs','MRs','lengthMicrons_smooth');
%     
%     pkmatun_delkaiBC = cat(1,pkmatun_delkaiBC,pkmatun);
% end


%% WT
%---- Peak amplitude of detrended data -----%
Amp = pkmatun_wt(:,5); % Peak amplitude of detrended data
CT = pkmatun_wt(:,4); % Circadian time
%SB = pkmatun_wt(:,9); % Cell length at birth

clear xs1 yperx1 lperx1

xs1 = unique(CT);
yperx1 = cell(length(xs1),1);
lperx1 = cell(length(xs1),1);

 for j = 1:length(CT)
        whereinxs = find(xs1 == CT(j));
        yperx1{whereinxs} = cat(1,yperx1{whereinxs},Amp(j));
        %lperx1{whereinxs} = cat(1,lperx1{whereinxs},SB(j));
 end

 clear meanyperx1 stdyperx1 serrmean1
 
for i = 1:length(yperx1)
    meanyperx1(i) = mean(yperx1{i});
    %meanlperx1(i) = mean(lperx1{i});
    stdyperx1(i) = std(yperx1{i});
    serrmean1(i) = stdyperx1(i)./sqrt(length(yperx1{i}));
end

clear xs2 meanyperx2 stdyperx2

xs2 = cat(1,xs1,xs1+15);
meanyperx2 = [meanyperx1,meanyperx1];
stdyperx2 = [stdyperx1,stdyperx1];

clear curve1 curve2 inBetween h 

figure
%hold on;
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
%hold on;plot(xs2,smooth(meanyperx2,9),'b-','linewidth',2)
set(gca,'XTick',7.5*(0:40));
grid on
nightmodeon15(s{22,2},7.5,1);
xlim([0 30])
ylim([0 500])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Pulse amplitude (a.u.)','FontWeight','bold')
%title('WT')
set(gca,'FontSize',14)
%saveas(gcf,['Fig5_4A'],'fig');saveas(gcf,['Fig5_4A'],'svg');

%% Plot Fig3D (no bining normalised)
meany = mean(meanyperx1);

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

%% binning 1.5 h normalised
meany = mean(meanyperx1);
meanyperx1_bin = [];
xs1_bin = [];

for i = 1:size(meanyperx1,2)/2
    meanyperx1_bin(i) = (meanyperx1(i*2-1)+meanyperx1(i*2))/2;
    xs1_bin(i) = (xs1(i*2-1)+xs1(i*2))/2;
end

xs2_bin = [xs1_bin,xs1_bin+15];
meanyperx2_bin = [meanyperx1_bin,meanyperx1_bin];


figure;
plot(xs2_bin, meanyperx2_bin./meany, '-','color',[0.8500 0.3250 0.0980],'LineWidth', 2);
%nightmodeon165(s{22,2},8.25,1);
xlim([0 30])
%ylim([0.5 2])
set(gca,'XTick',7.5*(0:40));
grid on
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Norm. pulse amplitude','FontWeight','bold')

%% WT
%---- Peak amplitude of detrended data -----%
Amp = pkmatun_wt(:,5); % Peak amplitude of detrended data
T = pkmatun_wt(:,2); % Time
%SB = pkmatun_wt(:,9); % Cell length at birth

clear xs1 yperx1 lperx1

xs1 = unique(T);
yperx1 = cell(length(xs1),1);
lperx1 = cell(length(xs1),1);

 for j = 1:length(T)
        whereinxs = find(xs1 == T(j));
        yperx1{whereinxs} = cat(1,yperx1{whereinxs},Amp(j));
        %lperx1{whereinxs} = cat(1,lperx1{whereinxs},SB(j));
 end

 clear meanyperx1 stdyperx1 serrmean1
 
for i = 1:length(yperx1)
    meanyperx1(i) = mean(yperx1{i});
    %meanlperx1(i) = mean(lperx1{i});
    stdyperx1(i) = std(yperx1{i});
    serrmean1(i) = stdyperx1(i)./sqrt(length(yperx1{i}));
end


clear curve1 curve2 inBetween h 

figure
%hold on;
curve1 = meanyperx1 + stdyperx1;
curve2 = meanyperx1 - stdyperx1;
x2 = [xs1', fliplr(xs1')];
inBetween = [curve1, fliplr(curve2)];
%h = fill(x2, inBetween,[0.8500, 0.3250, 0.0980],'LineStyle','none');
h = fill(x2, inBetween,'c','LineStyle','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',0.2)
hold on;
% plot(xs2, meanyperx2, '-','color',[0.8500, 0.3250, 0.0980],'LineWidth', 2);
plot(xs1, meanyperx1, '-','color','c','LineWidth', 2);
%hold on;plot(xs1,smooth(meanyperx1,9),'b-','linewidth',2)
hold on;
scatter(pkmatun_wt(:,2),pkmatun_wt(:,5),10,'k.')
set(gca,'XTick',8*(0:40));
grid on
nightmodeon16(schnitz2{2,2},8,1);
xlim([0 72])
ylim([0 500])
xlabel('Time in LL (h)','FontWeight','bold')
ylabel('Pulse amplitude (a.u.)','FontWeight','bold')
%title('WT')
set(gca,'FontSize',14)
%saveas(gcf,['Fig5_4A'],'fig');saveas(gcf,['Fig5_4A'],'svg');

%% WT
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
set(gca,'XTick',8.25*(0:40));
grid on
nightmodeon165(s{2,2},8.25,1);
xlim([0 49.5])
ylim([0.5 2.5])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Norm. pulse amplitude','FontWeight','bold')
%title('WT')
set(gca,'FontSize',18)
%saveas(gcf,['Fig5_4A'],'fig');saveas(gcf,['Fig5_4A'],'svg');

%%
%---- Peak width of detrended data -----%
Width = pkmatun_wt(:,6); % Peak width of detrended data
CT = pkmatun_wt(:,4); % Circadian time

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

xs2 = cat(1,xs1,xs1+16);
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
set(gca,'XTick',8*(0:40));
grid on
%nightmodeon(s{2,2},12,1);
xlim([0 32])
%ylim([1 7])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Pulse width (h)','FontWeight','bold')
%title('WT')
set(gca,'FontSize',18)

%% 

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
%%
figure; %% WT
%subplot(2,2,1);
%Y_wt = Y_wt-100;
for i = 1:size(Y_wt,2)
    p = patchline(T_wt(:,i),Y_wt(:,i),'linestyle','-','edgecolor','[0.1 0.1 0.1]','linewidth',0.5,'edgealpha',0.1);
    p.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold on
end
% hold on;
% [pk,lk] = findpeaks(rmmissing(Y_wt(:,1129)),rmmissing(T_wt(:,1129)),'MinPeakProminence',100);
% plot(T_wt(:,1129),Y_wt(:,1129),'m-',lk,pk,'mo','LineWidth',2, 'MarkerSize',8);
% hold on;
% [pk,lk] = findpeaks(rmmissing(Y_wt(:,1191)),rmmissing(T_wt(:,1191)),'MinPeakProminence',100);
% plot(T_wt(:,1191),Y_wt(:,1191),'r-',lk,pk,'ro','LineWidth',2, 'MarkerSize',8);
hold on;

% for i = 1:length(sall_wt)
%     sall_wt(i).MYs = sall_wt(i).MYs-100;
% end

%[xs, meanyperx, stdyperx, serrmean, ncells]= plotmeannew(sall_wt,'hrst','MYs',[],'c',1);
% [pk,lk] = findpeaks(meanyperx,xs,'MinPeakProminence',50);
% hold on
% plot(lk,pk,'mo','LineWidth',2, 'MarkerSize',6)%,'MarkerFaceColor','c');
set(gca,'XTick',7.5*(0:20));
grid on;
xlim([0 67.5])
ylabel('Dtr. P_{\itrpoD4}-EYFP (a.u.)','fontweight','bold')
xlabel('Time in LL (h)','fontweight','bold')
%title('Medium light (ccd = 10.5 h)')
nightmodeon15(s_ML{10,2},7.5,1);
ylim([-200 800])
set(gca,'FontSize',16)

%% WT
%---- Peak amplitude of detrended data over time -----%
Amp = pkmatun_wt(:,5); % Peak amplitude of detrended data
Ti = pkmatun_wt(:,2); % Time of peak

clear xs1 yperx1 lperx1

xs1 = unique(Ti);
yperx1 = cell(length(xs1),1);


 for j = 1:length(Ti)
        whereinxs = find(xs1 == Ti(j));
        yperx1{whereinxs} = cat(1,yperx1{whereinxs},Amp(j));
 end
 
clear meanyperx1 stdyperx1 serrmean1
 
for i = 1:length(yperx1)
    meanyperx1(i) = mean(yperx1{i});
    stdyperx1(i) = std(yperx1{i});
    serrmean1(i) = stdyperx1(i)./sqrt(length(yperx1{i}));
end


clear curve1 curve2 inBetween h 

hold on;
curve1 = (meanyperx1 + stdyperx1);
curve2 = (meanyperx1 - stdyperx1);
x2 = [xs1', fliplr(xs1')];
inBetween = [curve1, fliplr(curve2)];
%h = fill(x2, inBetween,[0.8500, 0.3250, 0.0980],'LineStyle','none');
h = fill(x2, inBetween,'c','LineStyle','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',0.2)
hold on;
% plot(xs2, meanyperx2, '-','color',[0.8500, 0.3250, 0.0980],'LineWidth', 2);
plot(xs1, meanyperx1, '-','color','c','LineWidth', 2);
set(gca,'XTick',7.5*(0:40));
grid on
%nightmodeon15(s_ML{20,2},7.5,1);
xlim([0 67.5])
ylim([-200 800])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Pulse amplitude (a.u.)','FontWeight','bold')
%title('WT')
set(gca,'FontSize',14)
%saveas(gcf,['Fig5_4A'],'fig');saveas(gcf,['Fig5_4A'],'svg');

%%
figure; %% WT
%subplot(2,2,1);
%Y_wt = Y_wt-100;
for i = 1:size(L_wt,2)
    p = patchline(T_wt(:,i),L_wt(:,i),'linestyle','-','edgecolor','[0.1 0.1 0.1]','linewidth',0.5,'edgealpha',0.1);
    hold on
end
% hold on;
% [pk,lk] = findpeaks(rmmissing(Y_wt(:,1129)),rmmissing(T_wt(:,1129)),'MinPeakProminence',100);
% plot(T_wt(:,1129),Y_wt(:,1129),'m-',lk,pk,'mo','LineWidth',2, 'MarkerSize',8);
% hold on;
% [pk,lk] = findpeaks(rmmissing(Y_wt(:,1191)),rmmissing(T_wt(:,1191)),'MinPeakProminence',100);
% plot(T_wt(:,1191),Y_wt(:,1191),'r-',lk,pk,'ro','LineWidth',2, 'MarkerSize',8);
hold on;

% for i = 1:length(sall_wt)
%     sall_wt(i).MYs = sall_wt(i).MYs-100;
% end

%[xs, meanyperx, stdyperx, serrmean, ncells]= plotmeannew(sall_wt,'hrst','MYs',[],'c',1);
% [pk,lk] = findpeaks(meanyperx,xs,'MinPeakProminence',50);
% hold on
% plot(lk,pk,'mo','LineWidth',2, 'MarkerSize',6)%,'MarkerFaceColor','c');
set(gca,'XTick',8*(0:20));
grid on;
%xlim([0 72])
ylabel('Cell length (\mum)','fontweight','bold')
xlabel('Time in LL (h)','fontweight','bold')
%title('Medium light (ccd = 10.5 h)')
%nightmodeon(s{10,2},12,1);
%ylim([200 1300])
set(gca,'FontSize',16)

%%
sall_R393C = [];

for k = 8:12
s1 = s{k,2};
sall_R393C = cat(2,sall_R393C,s1);
end

%%
sall_kaiR393C = [];

for k = 2:26
s1 = s{k,2};
sall_kaiR393C = cat(2,sall_kaiR393C,s1);
end
%%
[age, ctob,ctod,initial_size, final_size, size_diff, relsize_diff, elongrate, cellno] = schnitzphysprops2(sall_R393C,[],'hrst','lengthMicrons');

%%
figure;
histogram(ctod,0:1.5:24)
nightmodeon(s{9,2},12,1);
xlim([0 24])
ylim([0 70])
set(gca,'XTick',3*(0:40));
xlabel('Circadian time of division (h)','fontweight','bold')
ylabel('Frequency','fontweight','bold')

%%
figure;
histogram(age,0:1.5:36)
%nightmodeon(s{9,2},12,1);
xlim([0 36])
%ylim([0 70])
set(gca,'XTick',3*(0:40));
xlabel('Cell cycle duration (h)','fontweight','bold')
ylabel('Frequency','fontweight','bold')