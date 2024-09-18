function nightmodeon(schnitz,phase,rfp)

%tph = mod((min([schnitz.hrs]):0.01:max([schnitz.hrs]))-phase,25.5);
tph = mod((min([schnitz.hrs]):0.01:max([schnitz.hrs]))-phase,24);
isnight = zeros(1,length(tph));
for k = 1:length(tph)
    if tph(k)>=0 && tph(k)<12
        isnight(k)=1;
    else
        isnight(k)=0;
    end
end

hold on

if rfp == 1
    h = area(min([schnitz.hrs]):0.01:max([schnitz.hrs]),isnight*max([schnitz.MRs])*10,'FaceColor',[0.5,0.5,0.5],'LineStyle','none');
    h.FaceAlpha = 0.5;
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold on
    h1 = area(min([schnitz.hrs]):0.01:max([schnitz.hrs]),-isnight*max([schnitz.MRs])*10,'FaceColor',[0.5,0.5,0.5],'LineStyle','none');
    h1.FaceAlpha = 0.5;
    h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
else
    h = area(min([schnitz.hrs]):0.01:max([schnitz.hrs]),isnight*max([schnitz.MYs])*1.1,'FaceColor',[0.6,0.6,0.6],'LineStyle','none');
%     h = area(min([schnitz.hrs]):0.5:150,isnight*max([schnitz.MCs])*1.1,'FaceColor',[0.6,0.6,0.6],'LineStyle','none');
    h.FaceAlpha = 0.5;
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold on
    h1 = area(min([schnitz.hrs]):0.01:max([schnitz.hrs]),-isnight*max([schnitz.MYs])*1.1,'FaceColor',[0.5,0.5,0.5],'LineStyle','none');
    h1.FaceAlpha = 0.5;
    h1.Annotation.LegendInformation.IconDisplayStyle = 'off';

end

uistack(h,'bottom')
uistack(h1,'bottom')

hold off