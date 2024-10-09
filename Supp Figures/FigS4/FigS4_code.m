%% Fig.S4A (plot circadian time of birth vs birth length)
clear
load("FigS4AB_data.mat")

clear xs1 yperx1 xs2 meanyperx1 stdyperx1 meanyperx2 stdyperx2 s2

xs1 = unique(CTOB_delkaiBC);
yperx1 = cell(length(xs1),1);


 for j = 1:length(CTOB_delkaiBC)
        whereinxs = find(xs1 == CTOB_delkaiBC(j));
        yperx1{whereinxs} = cat(1,yperx1{whereinxs},LB_delkaiBC(j));
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

meanyperx1 = [];
stdyperx1 = [];
serrmean1 = [];

clear xs1 yperx1 xs2 meanyperx1 stdyperx1 meanyperx2 stdyperx2 s2

xs1 = unique(CTOB_delkaiBCdelrpoD4);
yperx1 = cell(length(xs1),1);


 for j = 1:length(CTOB_delkaiBCdelrpoD4)
        whereinxs = find(xs1 == CTOB_delkaiBCdelrpoD4(j));
        yperx1{whereinxs} = cat(1,yperx1{whereinxs},LB_delkaiBCdelrpoD4(j));
 end

for i = 1:length(yperx1)
    meanyperx1(i) = mean(yperx1{i});
    stdyperx1(i) = std(yperx1{i});
    serrmean1(i) = stdyperx1(i)./sqrt(length(yperx1{i}));
end

xs2 = cat(1,xs1,xs1+24);
meanyperx2 = [meanyperx1,meanyperx1];
stdyperx2 = [stdyperx1,stdyperx1];

hold on;
curve1 = meanyperx2 + stdyperx2;
curve2 = meanyperx2 - stdyperx2;
x2 = [xs2', fliplr(xs2')];
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween, 'y','LineStyle','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',0.2)
hold on;
plot(xs2, meanyperx2, '-','color','y', 'LineWidth', 2);


legend('\it\DeltakaiBC','\it\DeltakaiBC\DeltarpoD4')
set(gca,'XTick',12*(0:40));
grid on
xlim([0 48])
xlabel('Time at birth (h)','FontWeight','bold')
ylabel('Birth length (\mum)','FontWeight','bold')
set(gca,'FontSize',20)
nightmodeon(s{2,2},12,1);
ylim([2.5 4])


%% Fig.S4B (violin plot of birth length)
figure;
LB = [LB_WT;LB_delD4;LB_delkaiBC;LB_delkaiBCdelrpoD4];
groupIdx = [ones(size(LB_WT));2*ones(size(LB_delD4));3*ones(size(LB_delkaiBC));4*ones(size(LB_delkaiBCdelrpoD4))];
violin = violinplot(LB,groupIdx);
violin(1,1).ViolinColor = 'c';
violin(1,2).ViolinColor = 'b';
violin(1,3).ViolinColor = 'r';
violin(1,4).ViolinColor = 'y';
xticklabels({'WT','\it\DeltarpoD4','\it\DeltakaiBC','\it\DeltakaiBC\DeltarpoD4'});
ylabel('Birth length (\mum)','fontweight','bold')
ylim([2 5.5])
set(gca,'YTick',0.5*(0:40));
set(gca,'FontSize',14);


%% Fig.S4C (violin plot of cell cycle duration)
clear
load("FigS4CD_data.mat")

figure;
AGE = [AGE_WT;AGE_rpoD4OX];
groupIdx = [ones(size(AGE_WT));2*ones(size(AGE_rpoD4OX))];
violin = violinplot(AGE,groupIdx);
xticklabels({'WT','\itrpoD4OX'});
ylabel('Cell cycle duration (h)','fontweight','bold')
set(gca,'FontSize',16);

violin(1,1).ViolinColor = 'c';
violin(1,2).ViolinColor = 'm';
violin(1,1).ScatterPlot.SizeData = 10;
violin(1,2).ScatterPlot.SizeData = 10;


%% Fig.S4D (violin plot of elongation rate)
figure;
ER = [ER_WT;ER_rpoD4OX];
groupIdx = [ones(size(ER_WT));2*ones(size(ER_rpoD4OX))];
violin = violinplot(ER,groupIdx);
xticklabels({'WT','\itrpoD4OX'});
ylabel('Elongation rate (h^{-1})','fontweight','bold')
set(gca,'FontSize',16);

violin(1,1).ViolinColor = 'c';
violin(1,2).ViolinColor = 'm';
violin(1,1).ScatterPlot.SizeData = 10;
violin(1,2).ScatterPlot.SizeData = 10;
