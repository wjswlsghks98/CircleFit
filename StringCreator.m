function strs = StringCreator(varargin)
    % Create Strings for arc segment modes
    % Implemented by JinHwan Jeon, 2022

    %% INPUT
    % seg_type: Description of segment type: init, mid, last
    % cur_loc: Current index inside the segment of interest
    % n: Total number for small segments inside the segment of interest
    % prev_mode(not used for initial function call): Previous segment mode,
    % consists of "front","back","free","both", this wordings indicate
    % which side of the segment is fixed while performing Least Squares
    % fitting.
    %% OUTPUT
    % strs: String type array that consists of all possible segment modes
    %% OVERALL DESCRIPTION
    % All possible combination of segment modes are computed recursively
    % depending on the position of 'initial segment of interest'
    % After computing all the possible combination of segment modes, need
    % to find the order of arc spline approximation for every string in the
    % 'str' string array
    %% FUNCTION SCRIPT
    seg_type = varargin{1};
    cur_loc = varargin{2};
    n = varargin{3};
    switch seg_type
        case 'init'
            if n == 1
                strs = "Back";
            else
                % Initial Segment
                if cur_loc == 1
                    % Case 1: Start with "free"
                    strs1_ = StringCreator('init',cur_loc+1,n,"Free");
                    strs1 = [repmat("Free",size(strs1_,1),1) strs1_];
                    
                    % Case 2: Start with "back"
                    strs2_ = StringCreator('init',cur_loc+1,n,"Back");
                    strs2 = [repmat("Back",size(strs2_,1),1) strs2_];
                    
                    strs = [strs1; strs2];
                % Final Segment
                elseif cur_loc == n
                    prev_mode = varargin{4};
                    if strcmp(prev_mode,"Front") || strcmp(prev_mode, "Free")
                        strs = "Both";
                    else
                        strs = "Back";
                    end
                % Others
                else
                    prev_mode = varargin{4};
                    if strcmp(prev_mode,"Front") || strcmp(prev_mode,"Free")
                        strs1_ = StringCreator('init',cur_loc+1,n,"Front");
                        strs1 = [repmat("Front",size(strs1_,1),1) strs1_];
                        
                        strs2_ = StringCreator('init',cur_loc+1,n,"Both");
                        strs2 = [repmat("Both",size(strs2_,1),1) strs2_];
                        
                        strs = [strs1; strs2];
                        
                    elseif strcmp(prev_mode,"Back") || strcmp(prev_mode,"Both")
                        strs1_ = StringCreator('init',cur_loc+1,n,"Free");
                        strs1 = [repmat("Free",size(strs1_,1),1) strs1_];
                        
                        strs2_ = StringCreator('init',cur_loc+1,n,"Back");
                        strs2 = [repmat("Back",size(strs2_,1),1) strs2_];
                        
                        strs = [strs1; strs2];
                    end
                end
            end
        case 'last'
        case 'mid'
    end
end
