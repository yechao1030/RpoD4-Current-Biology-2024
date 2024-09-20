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

load('delkaiBC_ML_2475_DT_data.mat')

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