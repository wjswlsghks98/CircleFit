function order = FindSegOrder(strs)
    %  Finding Order of segment modes
    %% INPUT
    % strs: String type array that consists of all possible segment modes
    %% OUTPUT
    % order: Order of arc spline approximation for every string sequence
    % provided in "strs"
    %% OVERALL DESCRIPTION
    % All possible combination of segment modes are computed recursively
    % depending on the position of 'initial segment of interest'
    % After computing all the possible combination of segment modes, need
    % to find the order of arc spline approximation for every string in the
    % 'str' string array
    %% FUNCTION SCRIPT
    [m,n] = size(strs);
    order = zeros(m,n);
    for i=1:m
        j=1;
        accum_cnt = 0;
        accum_stack = [];
        while j <= n
            if strcmp(strs(i,j),"Front")
                accum_cnt = accum_cnt + 1;
                order(i,j) = accum_cnt;
            elseif strcmp(strs(i,j),"Back")
                accum_stack = [accum_stack j];
            elseif strcmp(strs(i,j),"Both")
                accum_stack = [accum_stack j];
            elseif strcmp(strs(i,j),"Free")
                accum_cnt =  accum_cnt + 1;
                order(i,j) = accum_cnt;
                % Pop accumulated stack
                L = length(accum_stack);

                for k=L:-1:1
                    accum_cnt = accum_cnt + 1;
                    order(i,accum_stack(k)) = accum_cnt;
                end
                accum_stack = [];
            end

            if ~isempty(accum_stack) && j == n
                % Pop accumulated stack
                L = length(accum_stack);

                for k=L:-1:1
                    accum_cnt = accum_cnt + 1;
                    order(i,accum_stack(k)) = accum_cnt;
                end
            end

            j = j+1;

        end
    end
end
