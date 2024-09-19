%% Very low light

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

%%
figure
CT = (0:0.75:23.75);
%cmp = parula(34);
cmp = parula(16);
cmp = cat(1,flipud(cmp),cmp);
for k = 1:size(T_wt,2)
    y = rmmissing(Y_wt(:,k));
    t = rmmissing(T_wt(:,k));
%     f = polyfit(t,y,2);
%     ypol = y-t.^2*f(1)-t*f(2)-f(3);
    [pks locs w p] = findpeaks(y,t,'MinPeakProminence',25);
    [pks0 idx0] = findpeaks(y,'MinPeakProminence',25);
    for n = 1:length(locs)
        ct = mod(locs(n),24);
        idx = find(CT==ct);
        %plot(t([max(idx0(n)-9,1):min(idx0(n)+9,length(y))]),y([max(idx0(n)-9,1):min(idx0(n)+9,length(y))]),'color',[cmp(idx,:)])
        plot(t([max(idx0(n)-9,1):min(idx0(n)+9,length(y))]),y([max(idx0(n)-9,1):min(idx0(n)+9,length(y))])./meany_VLL,'color',[cmp(idx,:)])% normalise by mean pulse height
        
        
        hold on
    end
end
set(gca,'XTick',12*(0:40));
grid on;
%nightmodeon(s_VLL{8,2},12.75,1);
xlim([0 120])
ylim([-0.5 2.5])
colormap(cmp)

colorbar('Ticks',[0:1/32:1],...
    'TickLabels',{'0','','','','3','','','','6','','','','9','','',...
    '','12','','','','15','','','','18','','','','21','','','','24'})

xlabel('Time in LL (h)','FontWeight','bold')
%ylabel('Detrended P_{\itrpoD4}-EYFP pulse (a.u.)','FontWeight','bold')
ylabel('Norm. detr. P_{\itrpoD4}-EYFP pulse','FontWeight','bold')
%title('Very low light (10 \muE m^{-2} s^{-1})')
title('10 \mumol m^{-2} s^{-1} (ccd = 22.5 h)')



%% Medium light WT
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
%%
figure
CT = (0:0.75:23.75);
%cmp = parula(32);
cmp = parula(16);
cmp = cat(1,flipud(cmp),cmp);
for k = 1:size(T_wt,2)
    y = rmmissing(Y_wt(:,k));
    t = rmmissing(T_wt(:,k));
%     f = polyfit(t,y,2);
%     ypol = y-t.^2*f(1)-t*f(2)-f(3);
%     [pks locs w p] = findpeaks(ypol,t,'MinPeakProminence',100);
%     [pks0 idx0] = findpeaks(ypol,'MinPeakProminence',100);
    [pks locs w p] = findpeaks(y,t,'MinPeakProminence',100);
    [pks0 idx0] = findpeaks(y,'MinPeakProminence',100);


    for n = 1:length(locs)
        ct = mod(locs(n),24);
        idx = find(CT==ct);
        %plot(t([max(idx0(n)-9,1):min(idx0(n)+9,length(y))]),y([max(idx0(n)-9,1):min(idx0(n)+9,length(y))]),'color',[cmp(idx,:)])
        plot(t([max(idx0(n)-9,1):min(idx0(n)+9,length(y))]),y([max(idx0(n)-9,1):min(idx0(n)+9,length(y))])./meany_ML,'color',[cmp(idx,:)])% normalise by mean pulse height
        
        
        hold on
    end
end
set(gca,'XTick',12*(0:40));
grid on;
%nightmodeon(s_ML{20,2},12,1);
xlim([0 96])
ylim([-1 2.5])
%ylim([-200 500])
colormap(cmp)
% colorbar('Ticks',[0:1/32:1],...
%     'TickLabels',{'','0','0.75','1.5','2.25','3','3.75','4.5','5.25','6','6.75','7.5','8.25','9','9.75','10.5',...
%     '11.25','12','12.75','13.5','14.25','15','15.75','16.5','17.25','18','18.75','19.5','20.25','21','21.75','22.5','23.25'})
colorbar('Ticks',[0:1/32:1],...
    'TickLabels',{'0','','','','3','','','','6','','','','9','','',...
    '','12','','','','15','','','','18','','','','21','','','','24'})
xlabel('Time in LL (h)','FontWeight','bold')
%ylabel('Detrended P_{\itrpoD4}-EYFP pulse (a.u.)','FontWeight','bold')
ylabel('Norm. detr. P_{\itrpoD4}-EYFP pulse','FontWeight','bold')
%title('Low light (18 \muE m^{-2} s^{-1})')
title('18 \mumol m^{-2} s^{-1} (ccd = 14.25 h)')



%% High light
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
%%
figure
CT = (0:0.75:23.75);
%cmp = parula(32);
cmp = parula(16);
cmp = cat(1,flipud(cmp),cmp);
for k = 1:size(T_wt,2)
    y = rmmissing(Y_wt(:,k));
    t = rmmissing(T_wt(:,k));
%     f = polyfit(t,y,2);
%     ypol = y-t.^2*f(1)-t*f(2)-f(3);
    [pks locs w p] = findpeaks(y,t,'MinPeakProminence',100);
    [pks0 idx0] = findpeaks(y,'MinPeakProminence',100);
    for n = 1:length(locs)
        ct = mod(locs(n),24);
        idx = find(CT==ct);
        %plot(t([max(idx0(n)-9,1):min(idx0(n)+9,length(y))]),y([max(idx0(n)-9,1):min(idx0(n)+9,length(y))]),'color',[cmp(idx,:)])
        plot(t([max(idx0(n)-9,1):min(idx0(n)+9,length(y))]),y([max(idx0(n)-9,1):min(idx0(n)+9,length(y))])./meany_HL,'color',[cmp(idx,:)])%normalised by mean pulse height
        
        hold on
    end
end
set(gca,'XTick',12*(0:40));
grid on;
%nightmodeon(s_VLL{8,2},12,1);
xlim([0 72])
ylim([-1 2.5])
colormap(cmp)
% colorbar('Ticks',[0:1/32:1],...
%     'TickLabels',{'','0','0.75','1.5','2.25','3','3.75','4.5','5.25','6','6.75','7.5','8.25','9','9.75','10.5',...
%     '11.25','12','12.75','13.5','14.25','15','15.75','16.5','17.25','18','18.75','19.5','20.25','21','21.75','22.5','23.25'})
colorbar('Ticks',[0:1/32:1],...
    'TickLabels',{'0','','','','3','','','','6','','','','9','','',...
    '','12','','','','15','','','','18','','','','21','','','','24'})
xlabel('Time in LL (h)','FontWeight','bold')
ylabel('Norm. detr. P_{\itrpoD4}-EYFP pulse','FontWeight','bold')
%title('Medium light (35 \muE m^{-2} s^{-1})')
title('25 \mumol m^{-2} s^{-1} (ccd = 10.5 h)')
