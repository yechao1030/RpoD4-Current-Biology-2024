function [pkmatun] = findpeaksmethod_rpoD4(schnitz,time,ylabel,length);

[Y salt] = schnitzTables(schnitz,ylabel);
[T salt] = schnitzTables(schnitz,time);
[L salt] = schnitzTables(schnitz,length);


cutoff = 0;

Pks = []; % Peak amplitude
Locs = []; % Location of the peak
CTLocs = []; % Circadian time of locations
W = []; % Width of the peak
P = []; % Prominence of the peak
Schnitz = []; % lineage
YFP = []; % amplitude before detrending
CL = []; % cell length at YFP peak

    for k = 1:size(T,2)
        t = T(:,k);
        y = Y(:,k);
        l = L(:,k);
        schnitzes = salt(:,k);
        y = y(t>=cutoff);
        s = salt(t>=cutoff);
        t = t(t>=cutoff);
        %f = polyfit(t,y,2); % detrend data by fitting a polynormial curve for the data
%        ypol = y-t.^2*f(1)-t*f(2)-f(3);
        if size(t,1)>=3
%         [pks locs w p] = findpeaks(ypol,t,'MinPeakProminence',100);
%         [pks0 idx0] = findpeaks(ypol,'MinPeakProminence',100);
        [pks locs w p] = findpeaks(y,t,'MinPeakProminence',100);
        [pks0 idx0] = findpeaks(y,'MinPeakProminence',100);
        
        ctlocs = mod(locs,24);
        schnitz = schnitzes(idx0);
        yfp = y(idx0);
        cl = l(idx0);

        YFP = cat(1,YFP,yfp);
        Locs = cat(1,Locs,locs);
        Schnitz = cat(1,Schnitz,schnitz);
        CTLocs = cat(1,CTLocs,ctlocs);
        Pks = cat(1,Pks,pks);
        W = cat(1,W,w);
        P = cat(1,P,p);
        CL = cat(1,CL,cl);
  
        end
    end

    pkmat = [YFP,Locs,Schnitz,CTLocs,Pks,W,P,CL];
    [a b c] = unique(pkmat(:,1:3),'rows','stable'); % find unique peaks
    pkmatun = pkmat(b,:);