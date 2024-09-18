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

%%
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

%%
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

%% Plot Figure 3A panels

load Fig3A_data.mat

figure
%subplot(3,1,1)
%for i = 1:size(Y_VLL,1)
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


figure
%subplot(3,1,2)
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


figure
%subplot(3,1,3)
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

%%
figure;
for i = 1:size(Y_HL,2)
    subplot(2,1,1)
y = rmmissing(Y_HL(:,i));
x = rmmissing(T_HL(:,i));
[pk,lk] = findpeaks(y,x,'MinPeakProminence',100);
plot(x,y/meany_HL,'-','Color',[0.85,0.33,0.10],'LineWidth',2);
hold on
plot(lk,pk./meany_HL,'o','Color',[0.85,0.33,0.10],'LineWidth',2, 'MarkerSize',8);
xlim([0 72])
ylim([-0.5 2])
set(gca,'XTick',12*(0:40));
grid on
xlabel('Time in LL (h)','fontweight','bold')
title('25 \muE m^{-2} s^{-1}',num2str(i))
hold off

subplot(2,1,2)
y1 = rmmissing(L_HL(:,i));
findpeaks(y1,x)
xlim([0 72])

set(gca,'XTick',12*(0:40));
grid on
pause
end


%%
figure;
subplot(3,1,1)
y = rmmissing(L_VLL(:,100));
x = rmmissing(T_VLL(:,100));
[pk,lk] = findpeaks(y,x,'MinPeakProminence',0.1);
plot(x,y,'c-',lk,pk,'mo','LineWidth',2, 'MarkerSize',12);
xlim([0 72])
set(gca,'XTick',12*(0:40));
grid on

subplot(3,1,2)
y = rmmissing(L_ML(:,300));
x = rmmissing(T_ML(:,300));
[pk,lk] = findpeaks(y,x,'MinPeakProminence',0.1);
plot(x,y,'c-',lk,pk,'mo','LineWidth',2, 'MarkerSize',12);
xlim([0 72])
set(gca,'XTick',12*(0:40));
grid on

subplot(3,1,3)
y = rmmissing(L_HL(:,600));
x = rmmissing(T_HL(:,600));
[pk,lk] = findpeaks(y,x,'MinPeakProminence',0.1);
plot(x,y,'c-',lk,pk,'mo','LineWidth',2, 'MarkerSize',12);
xlim([0 72])
set(gca,'XTick',12*(0:40));
grid on
