function [data, schnitz_nums] = schnitzTables(s, channel)
% Function produces two matrices as an output. Each row corresponds to a
% frame, each column to a schnitz that has no further children. Data is
% backtracked up the tree so all parent data is included in a column.

% The first matrix contains the value of the selected channel. If the
% channel data does have the same length as the frame numbers for that
% schnitz, it is padded with NaNs (e.g. with the output of
% promoterActivity)

% The second matrix contains the source schnitz for the data in the first.
% This is useful to keep track of where the data is from and perhaps weight
% repeated occurences of parents.

    end_frame = max([s.frames]);
    start_frame = min([s.frames]);
    childless = ~([s.D] > 0 | [s.E] > 0);
    data = nan(end_frame, nnz(childless));
    schnitz_nums = nan(end_frame, nnz(childless));
    schnitzes = find(childless);
    
    mismatches = [];

    for i = 1:length(schnitzes)
        j = schnitzes(i);
        while j > 0
            schnitz_data = s(j).(channel);
            % Deal with mismatching frames/data lengths from
            % promoterActivity code
            if length(s(j).frames) ~= length(schnitz_data)
                if ~any(mismatches == j)
                    fprintf('WARNING: Frame & channel mismatch on schnitz %i\n', j)
                    fprintf('first found processing schnitz %i\n', schnitzes(i))
                    mismatches = [mismatches j];
                end
                len_diff = length(s(j).frames) - length(s(j).(channel));
                schnitz_data = padarray(schnitz_data, [0 len_diff], NaN, 'pre');
            end
            data(s(j).frames, i) = schnitz_data;
            schnitz_nums(s(j).frames, i) = j;
            j = s(j).P;
        end
    end
    
    % Delete empty frames. Still necessary to have them so indexing works
    % with frame numbers.
    data(1:(start_frame-1),:) = [];
    schnitz_nums(1:(start_frame-1),:) = [];
end