%% Fig.1C (WT)
load('2021-07-20_cyano-03_6.mat')

figure;

subplot(2,1,1);
plot(T,Y,'k-');
hold on;
[pk,lk] = findpeaks((Y(:,9)),T(:,9),'MinPeakProminence',100);
plot(T(:,9),Y(:,9),'c-',lk,pk,'mo','LineWidth',2, 'MarkerSize',12);
nightmodeon(s,12,1);
xlabel('')
ylabel('P_{rpoD4}-EYFP (a.u.)','fontweight','bold')
ylim([400 1200])
xlim([24 96])
set(gca,'XTick',12*(0:20));
set(gca,'FontSize',12);
grid on;

subplot(2,1,2);
plot(T,L_smooth,'k-');
hold on
TF = islocalmin(L_smooth(:,9));
x = T(:,9);y = L_smooth(:,9);
plot(x,y,'c-',x(TF),y(TF),'mo','LineWidth',2, 'MarkerSize',12);
nightmodeon(s,12,1);
xlabel('Time in LL (h)','fontweight','bold')
ylabel('Cell length (\mum)','fontweight','bold')
ylim([2 10])
xlim([24 96])
set(gca,'XTick',12*(0:20));
set(gca,'FontSize',12);
grid on;




%% Fig.1D_WT PrpoD4-EYFP vs ToD
clear

load('WT_ML_24_DT_data.mat')

%%%% Find individual traces of cells with complete cell cycles %%%%
T_wt = [];
T_wt_all = [];
Y_wt = [];
Y_wt_all = [];


for j = 2:length(s_ML)
    for i = 1:length(s_ML{j,2})
        if s_ML{j,2}(i).E ~= 0 && s_ML{j,2}(i).P ~= 0 % find cells with complete cell cycle
            A = s_ML{j,2}(i).E; % daughter1
            B = s_ML{j,2}(i).D; % daughter2
            t1 = [s_ML{j,2}(i).hrs-max(s_ML{j,2}(i).hrs)-0.75,s_ML{j,2}(A).hrs-min(s_ML{j,2}(A).hrs)]; % adjust the time according to time of division
            y1 = [s_ML{j,2}(i).MYs_det,s_ML{j,2}(A).MYs_det];
            
            
            t2 = [s_ML{j,2}(i).hrs-max(s_ML{j,2}(i).hrs)-0.75,s_ML{j,2}(B).hrs-min(s_ML{j,2}(B).hrs)];
            y2 = [s_ML{j,2}(i).MYs_det,s_ML{j,2}(B).MYs_det];
                   
           
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


%%%% Calculate mean trace from cells with complete cell cycles %%%%
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


%%%% Make the plot %%%%`
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
ylabel('P_{rpoD4}-EYFP (a.u.)','FontWeight','bold')
grid on;
xline(0,'k--','LineWidth',2)
xlim([-9,9]);
ylim([-200,600]);
set(gca,'XTick',3*(-40:40));
set(gca,'FontSize',14)


%% Fig.1E_delkaiBC PrpoD4-EYFP vs ToD
clear 

load('Fig1G_data.mat')

T_delkaiBC = [];
T_delkaiBC_all = [];
Y_delkaiBC = [];
Y_delkaiBC_all = [];


for j = 2:length(s_delkaiBC)
    for i = 1:length(s_delkaiBC{j,2})
        if s_delkaiBC{j,2}(i).E ~= 0 && s_delkaiBC{j,2}(i).P ~= 0
            A = s_delkaiBC{j,2}(i).E;
            B = s_delkaiBC{j,2}(i).D;
            t1 = [s_delkaiBC{j,2}(i).hrs-max(s_delkaiBC{j,2}(i).hrs)-0.75,s_delkaiBC{j,2}(A).hrs-min(s_delkaiBC{j,2}(A).hrs)];
            y1 = [s_delkaiBC{j,2}(i).MYs_det,s_delkaiBC{j,2}(A).MYs_det];
            
            t2 = [s_delkaiBC{j,2}(i).hrs-max(s_delkaiBC{j,2}(i).hrs)-0.75,s_delkaiBC{j,2}(B).hrs-min(s_delkaiBC{j,2}(B).hrs)];
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
ylabel('P_{rpoD4}-EYFP (a.u.)','FontWeight','bold')
grid on;
xline(0,'k--','LineWidth',2)
xlim([-9,9]);
ylim([-200,600]);
set(gca,'XTick',3*(-40:40));
set(gca,'FontSize',14)


%% Fig.1F WT
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

s = s_ML;

pkmatun_wt = [];

for i = 2:20
    [pkmatun] = findpeaksmethod_rpoD4(s{i,2},'hrst','MYs_det','lengthMicrons_smooth',24,0,100);
    
    pkmatun_wt = cat(1,pkmatun_wt,pkmatun);
end



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

%% Fig.1G delkaiBC
clear

load('Fig1G_data.mat')

for i = 2:19
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

s = s_delkaiBC;

pkmatun_delkaiBC = [];


for i = 2:19
    [pkmatun] = findpeaksmethod_rpoD4(s{i,2},'hrst','MYs_det','lengthMicrons_smooth',24,0,100);
    
    pkmatun_delkaiBC = cat(1,pkmatun_delkaiBC,pkmatun);
end


%---- Peak amplitude of detrended data -----%
Amp = pkmatun_delkaiBC(:,5); % Peak amplitude of detrended data
CT = pkmatun_delkaiBC(:,4); % Circadian time

clear xs1 yperx1 xs2 meanyperx2 stdyperx2 x2 inBetween
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
h = fill(x2, inBetween,'r','LineStyle','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',0.2)
hold on;
plot(xs2, meanyperx2, '-','color','r','LineWidth', 2);
set(gca,'XTick',12*(0:40));
grid on
nightmodeon(s_delkaiBC{2,2},12,1);
xlim([0 48])
ylim([100 450])
xlabel('Time of day (h)','FontWeight','bold')
ylabel('Pulse amplitude (a.u.)','FontWeight','bold')
title('\it\DeltakaiBC')
set(gca,'FontSize',18)