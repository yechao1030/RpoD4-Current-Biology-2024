function [pkmatun] = findpeaksmethod_rpoD4(schnitz,time,ylabel,length,period,cutoff_start,cutoff_end,prominence)

% function to find unique peaks of rpoD4 expression time traces
% "schnitz, time, ylabel & lenghth" are all Schnitzcell variables
% "period" is the clock period, the period for WT clock is 24, and the
% period for KaiR393C mutant is 15.
% "cutoff_start" & "cutoff_end" defines the time frame in which the analysis is
% carried out
% "prominence" is the threshold for filtering a genuine expression peak from
% background noise. Please see Methods for which prominence threshold is
% used for data acquired under different conditions

[Y salt] = schnitzTables(schnitz,ylabel);
[T salt] = schnitzTables(schnitz,time);
[L salt] = schnitzTables(schnitz,length);


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
        y = y(t>=cutoff_start & t<=cutoff_end);
        s = salt(t>=cutoff_start & t<=cutoff_end);
        t = t(t>=cutoff_start & t<=cutoff_end);

        if size(t,1)>=3
        [pks locs w p] = findpeaks(y,t,'MinPeakProminence',prominence);
        [pks0 idx0] = findpeaks(y,'MinPeakProminence',prominence);
        
        ctlocs = mod(locs,period);
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