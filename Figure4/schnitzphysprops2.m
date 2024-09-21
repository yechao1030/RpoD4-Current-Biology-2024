function [age, ctob,ctod,yob,initial_size, final_size, size_diff, relsize_diff, elongrate, cellno] = schnitzphysprops2(schnitz,whichschnitz,timefield,sizefield,fluorfield)
% function to compute distributions of single cell physiological properties (essentially cell cycle related properties) from the schnitz data structure
% The user can specify which schnitzes should be evaluated

if isempty(whichschnitz)
    % if argument is specified as [], then include all
    whichschnitz = 1:length(schnitz);
end

% the properties are cell cycle duration (age), size at division, size
% at start of cell cycle, absolute and relative added size during cell cycle
% the user can specify which schnitz field to consider as "size" and "time".

age = [];
ctob = [];
ctod = [];

yob = [];

initial_size = [];
final_size = [];
size_diff = [];
relsize_diff = [];
elongrate = [];
cellno = [];


for k = 1:length(whichschnitz)
    j = whichschnitz(k);
    
    if schnitz(j).D > 0 && schnitz(j).(timefield)(1) > 0 % if they have daughters and it's not the first cell cycle (which I'll exclude)
        size = schnitz(j).(sizefield);
        %final_size = cat(2,final_size,size(end));
        
        if schnitz(j).P > 0 % and if they complete a cell cycle with identified daughters
            agepre = schnitz(j).(timefield)(end)-schnitz(j).(timefield)(1);
            if agepre > 0
                age = cat(2,age,schnitz(j).(timefield)(end)-schnitz(j).(timefield)(1));
                ctob = cat(2,ctob,mod(schnitz(j).(timefield)(1),24));
                ctod = cat(2,ctod,mod(schnitz(j).(timefield)(end),24));
                yob = cat(2,yob,schnitz(j).(fluorfield)(1));
                size_diff = cat(2,size_diff,size(end)-size(1));
                relsize_diff = cat(2,relsize_diff,100*(size(end)-size(1))/size(1));
                initial_size = cat(2,initial_size,size(1));
                final_size = cat(2,final_size,size(end));
                % Elongation rate
                f1 = fit(schnitz(j).(timefield)',schnitz(j).(sizefield)','exp1');
                % f1 is fit structure with elements a and b, such that
                % y = a*exp(b*t). Therefore:
                elongrate = cat(2,elongrate,f1.b);
                cellno = cat(2,cellno,j);
                
            end
        end
    end
    
    if schnitz(j).P > 0 % if they have parents, regardless of whether they complete a full cycle
        size = schnitz(j).(sizefield);
        %initial_size = cat(2,initial_size,size(1));
    end
    
end
end

 
            