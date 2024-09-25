%% Parameters and equations used in Figure 2

clear 

tfreq = 0.01;
t = 0:tfreq:100;


freq1 = 12; % Change accordingly for Figure 2D
g1 = 50; % gc
a1 = 1; % ac
w1 = 2*pi/freq1; %  ωc average cell division frequency.
n1 = 4; % nc: Hill coefficient of cell division-driven signal.
k1 = 10; % kc
y1 = a1*(1 + cos(w1*t)); % Eq S4
g2 = 1; % gm
a2 = 1; % am
w2 = 2*pi/24; %  ωm Circadian frequency.
n2 = 1; % nm
k2 = 1; % km
y2 = a2*(1 + cos(w2*(t+12))); % Eq S5

u1 = (y1.^n1*g1)./(k1^n1 + y1.^n1);
u2 = 1+(y2.^n2*g2)./(1 + y2.^n2);

s = u1.*u2; % Eq S6

[pks, locs, w, prom] = findpeaks(s,t);

%% Fig2B_1
figure;
subplot(2,1,1)
plot(t,u1,'k-','linewidth',2)
set(gca,'XTick',0:6:100),grid on
xlabel('Time (h)','FontWeight','bold'), ylabel({'(Carrier)';'Cell cycle pulsing (a.u.)'},'FontWeight','bold') % legend('cell cycle pulsing','Clock')

subplot(2,1,2)
plot(t,u2,'r-','linewidth',2)
xlabel('Time (h)','FontWeight','bold'),ylabel({'(Modulator)';'Clock (a.u.)'},'FontWeight','bold')
set(gca,'XTick',0:6:100),grid on

%% Fig2B_2
figure;
findpeaks(s,t)
hold on
plot(t,s,'b-','linewidth',2)
xlabel('Time (h)','FontWeight','bold')
ylabel('RpoD4 (a.u.)','FontWeight','bold')
set(gca,'XTick',0:6:100),grid on

%% Fig2C
figure

% Loop through different starting conditions to produce color-coded peaks
cmp = parula(12);
cmp = cat(1,flipud(cmp),cmp);
for j = 0:freq1
    
    y1 = a1*(1 + cos(w1*(t+j)));
    u1 = (y1.^n1*g1)./(k1^n1 + y1.^n1);
    s = u1.*u2;
    [pks, locs, w, prom] = findpeaks(s,t);
    for k = 1:length(pks)
        locs_circ = round(mod(locs(k),24));
        twin = 1/2*freq1;
        tb = max(0,locs(k) - twin); tb = round(tb); indb = find(t == tb);
        tf = min(t(end),locs(k) + twin); tf = round(tf); indf = find(t == tf);
        plot(t(indb:indf), s(indb:indf), 'Color',[cmp(locs_circ+1,:)]),hold on
        
    end
end
xlabel('Time (h)','FontWeight','bold'),ylabel('RpoD4 (a.u.)','FontWeight','bold')
set(gca,'XTick',0:6:100),grid on,
hold off
colormap(cmp)
colorbar('Ticks',0:3:24),caxis([0 24])

%% Fig2D

x = 0:0.01:100;
ref1 = a1*(1 + cos(w1*0));
ref1 = (ref1.^n1*g1)./(k1^n1 + ref1.^n1);
for j = 1:length(x)
    ref2 = a2*(1 + cos(w2*(24+12+x(j))));
    ref2 = 1+(ref2.^n2*g2)./(1 + ref2.^n2);
    pam(j) = ref1.*ref2;
end

figure;
findpeaks(s,t)
hold on
plot(t,s,'b-','linewidth',2)
plot(x,pam,'c-','linewidth',2)
set(gca,'XTick',0:6:100)
xlim([0 100])
grid on
xlabel('RpoD4 (a.u.)','FontWeight','bold')
ylabel('Pulse amplitude (a.u.)','FontWeight','bold')

%% Fig2F
x = 0:0.01:100;
ref1 = a1*(1 + cos(w1*0));
ref1 = (ref1.^n1*g1)./(k1^n1 + ref1.^n1);
for j = 1:length(x)
    ref2 = a2*(1 + cos(w2*(24+12+x(j))));
    ref2 = 1+(ref2.^n2*g2)./(1 + ref2.^n2);
    pam(j) = ref1.*ref2;
end

figure
plot(x,pam,'c-','linewidth',4)
set(gca,'XTick',0:6:48)
xlim([0 48])
grid on
xlabel('Time (h)','FontWeight','bold')
ylabel('Pulse amplitude (a.u.)','FontWeight','bold')
