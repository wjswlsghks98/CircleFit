classdef CurveFitV2 < handle
    % CurveFitting for optimized point-wise dataset
    % Optimization done in Batch Manner(Not Incremental)
    % Implemented by JinHwan Jeon, 2022
    properties
        plot_num = 1 % Plot Number
        Optimizer
        opt = struct()
        
    end
    methods
        %% Constructor
        function obj = CurveFitV2(Optimizer)
            obj.Optimizer = Optimizer;
        end
        
        %% Part1 
        %% Part 1 Phase 1 (Modified): Find straight line segment candidates
        function obj = optimizePt1Ph1(obj)
            % Instead of performing Dynamic programming, naively search for
            % straight line segments
            LP_l = obj.Optimizer.opt.reordered_lml_pc(1:2,:);
            LP_r = obj.Optimizer.opt.reordered_lmr_pc(1:2,:);
            
            W_l = obj.Optimizer.opt.w_l;
            W_r = obj.Optimizer.opt.w_r;

            base_thres = 1000;
            lb = 1;
            
            plot_flag = false;
            
            obj.opt.line_seg_length = zeros(1,size(LP_l,2)-base_thres);
            obj.opt.line_intvs = [];

            while lb <= size(LP_l,2) - base_thres
                
                valid = true;
                ub = lb + base_thres;

                while valid
                    X_l = LP_l(1,lb:ub); Y_l = LP_l(2,lb:ub);
                    X_r = LP_r(1,lb:ub); Y_r = LP_r(2,lb:ub);
                    w_l = W_l(lb:ub); w_r = W_r(lb:ub);
                    [~,err] = LineFitV2(X_l, X_r, Y_l, Y_r, w_l, w_r, plot_flag);
                    
                    err_list = [err.full_l err.full_r];
                    if length(find(err_list > 5e-2)) >= 2
                        valid = false;
                    end
                    ub = ub + 1;
                end

                if (ub-1) - lb == base_thres
                    obj.opt.line_seg_length(lb) = 0;
                    if lb ~= 1 && obj.opt.line_seg_length(lb-1) ~= 0
                        obj.opt.line_intvs(end,2) = lb-1;
                    end
                else
                    obj.opt.line_seg_length(lb) = (ub-1) - lb;
                    if obj.opt.line_seg_length(lb-1) == 0
                        obj.opt.line_intvs = [obj.opt.line_intvs; [lb 0]];
                    end
                    disp(['Lower Bound Idx: ',num2str(lb), ' Maximum Segment Length: ',num2str((ub-1) - lb)])
                end
                
                lb = lb + 1;
            end

            obj.opt.line_seg_length = [obj.opt.line_seg_length zeros(1,base_thres)];
        end

        %% Part 1 Phase 1 (Modified): Extract straight line segments
        function obj = optimizePt1Ph2(obj)
            LP_l = obj.Optimizer.opt.reordered_lml_pc(1:2,:);
            LP_r = obj.Optimizer.opt.reordered_lmr_pc(1:2,:);
            
            W_l = obj.Optimizer.opt.w_l;
            W_r = obj.Optimizer.opt.w_r;
            plot_flag = false;
            seg_length = obj.opt.line_seg_length;

            obj.opt.line_segments = {};
            % From the maximum length segment, downward search
            [max_L,lb] = max(seg_length);
            obj.opt.line_segments = {};
            cnt = 1;
            arc_spacing = 100;
            min_line_length = 500;

            while max_L > 0
                
                ub = lb + seg_length(lb);

                if ub > length(seg_length)
                    seg_length(lb:end) = -1;
                else
                    if seg_length(ub) == -1
                        while seg_length(ub) == -1
                            ub = ub - 1;
                        end
                    end
                end
                
                
                if ub - lb < min_line_length
                    % If segment is too short, ignore
                    seg_length(lb:ub) = -1;
                    [max_L,lb] = max(seg_length);
                else
                    X_l = LP_l(1,lb:ub); Y_l = LP_l(2,lb:ub);
                    X_r = LP_r(1,lb:ub); Y_r = LP_r(2,lb:ub);
                    w_l = W_l(lb:ub); w_r = W_r(lb:ub);
                    [res,err] = LineFitV2(X_l, X_r, Y_l, Y_r, w_l, w_r, plot_flag);
    
                    res.status = 'line';
                    res.bnds = [lb ub];
                    obj.opt.line_segments = [obj.opt.line_segments {res}];
                    disp(['Line Segment ', num2str(cnt),': Idx ',num2str(lb),' ~ ', num2str(ub)])
                    disp(['Max Error: ', num2str(max([err.max_l err.max_r])), 'm RMSE Error: ', num2str(err.rmse), 'm'])
                    disp('----------------------------------------------')
                    
                    if ub <= length(seg_length)
                        seg_length(lb-arc_spacing:ub) = -1;
                    end
                    
                    if ub + arc_spacing <= length(seg_length)
                        seg_length(ub:ub+arc_spacing) = -1;
                    else
                        seg_length(ub:end) = -1;
                    end
                    
                    
                    [max_L,lb] = max(seg_length);
                    cnt = cnt + 1;
                end

            end

            % Need to reorder segments and create interval information
            n = length(obj.opt.line_segments);
            lbs = zeros(n,1);
            intvs = zeros(n,2);

            for i=1:n
                seg = obj.opt.line_segments{i};
                lbs(i) = seg.bnds(1);
                intvs(i,:) = seg.bnds;
            end
            
            [~,I] = sort(lbs);
            intvs = intvs(I,:);
            obj.opt.line_intvs = intvs;
            
            cpyd = obj.opt.line_segments;

            for i=1:n
                idx = I(i);
                obj.opt.line_segments{i} = cpyd{idx};
            end

        end

        %% Part2: Optimize
        %% Phase 2: Fill in remaining parts with circular arcs
        function obj = optimizePt2(obj)
            
            LP_l = obj.Optimizer.opt.reordered_lml_pc(1:2,:);
            LP_r = obj.Optimizer.opt.reordered_lmr_pc(1:2,:);
            
            W_l = obj.Optimizer.opt.w_l;
            W_r = obj.Optimizer.opt.w_r;

            disp('-Phase 2: Filling in Remaining Parts with Circular Arcs')
            
            n = size(obj.opt.line_intvs,1);
            obj.opt.arc_intvs = zeros(n+1,2);
            obj.opt.arc_intvs(1,1) = 1;
            obj.opt.arc_intvs(1,2) = obj.opt.line_intvs(1,1)-1;

            for i=1:n-1
                obj.opt.arc_intvs(i+1,1) = obj.opt.line_intvs(i,2)+1;
                obj.opt.arc_intvs(i+1,2) = obj.opt.line_intvs(i+1,1)-1;
            end

            obj.opt.arc_intvs(end,1) = obj.opt.line_intvs(end,2)+1;
            obj.opt.arc_intvs(end,2) = size(obj.Optimizer.opt.reordered_lml_pc,2);
            
            %% One-way optimization for first and last segments
            disp('-Arc Approximation for first segment-')
            
            seg = obj.opt.line_segments{1};
            
            obj.opt.arc_segments = {};
            err_thres = 0.15;
            err_max = inf;
            num_seg = 1;
            
            while err_max > err_thres
                disp('===============================================')
                disp(['Number of arc segments used for optimization: ',num2str(num_seg)])
                seg_strs = StringCreator('init',1,num_seg);
                seg_order = FindSegOrder(seg_strs);
                [~,n] = size(seg_order);                

                bnds = obj.opt.arc_intvs(1,:);
                if n == 1
                    lb_ = 1; ub_ = seg.bnds(1)-1;
                    X_l = LP_l(1,lb_:ub_); X_r = LP_r(1,lb_:ub_);
                    Y_l = LP_l(2,lb_:ub_); Y_r = LP_r(2,lb_:ub_);
                    w_l = W_l(lb_:ub_); w_r = W_r(lb_:ub_);
                    fixed_points = [seg.init_pointL seg.init_pointR];
                    [~,err] = CircleFitV3(X_l,X_r,Y_l,Y_r,w_l,w_r,...
                                          false,"Back",fixed_points,true);
                    err_max = err.max;
                    opt_bnds = [lb_ ub_];
                else
                    obj.opt.opt_ = FindBestIntvs(LP_l,LP_r,W_l,W_r,bnds,seg_strs,seg_order,{seg},'init');
                    opt_bnds = sort([1,seg.bnds(1)-1, obj.opt.opt_.bnds]);
                    err_max = max(obj.opt.opt_.err_maxs);
                end
                disp('Current optimized boundaries:')
                disp(opt_bnds)
                disp(['Maximum error for current segmentation: ',num2str(err_max),'m'])
                num_seg = num_seg + 1;
                
            end

            %% Two-way optimization for segments in between
            disp('-Arc Approximation for remaining segment-')

            
        end

        %% Phase 3: Shift Segments to match lateral position (Optimization)
        function obj = optimizePh3(obj)
        end
        
        %% Plot Results
        function plotRes(obj)
            line_segs = obj.opt.line_segments;
%             arc_segs = obj.opt.arc_segments;
            
            LP_l = obj.Optimizer.opt.reordered_lml_pc(1:2,:);
            LP_r = obj.Optimizer.opt.reordered_lmr_pc(1:2,:);

            % Plot Line Segments 
            figure(2); hold on; axis equal; grid on;
            plot(LP_l(1,:),LP_l(2,:),'r.');
            plot(LP_r(1,:),LP_r(2,:),'g.');

            n = length(line_segs);
            for i=1:n
                seg = line_segs{i};
                left_lane = [seg.init_pointL seg.last_pointL];
                right_lane = [seg.init_pointR seg.last_pointR];
              
                plot(left_lane(1,:),left_lane(2,:),'k--')
                plot(right_lane(1,:),right_lane(2,:),'k--')
                plot(left_lane(1,:),left_lane(2,:),'bp')
                plot(right_lane(1,:),right_lane(2,:),'bp')
            end

            % Plot Arc Segments
%             n = length(arc_segs);
%             for i=1:n
%                 seg = arc_segs{i};
%                 ths = linspace(seg.th_init, seg.th_last, 1e4);
%                 
%                 R = seg.R; delL = seg.delL;
%                 x = seg.x; y = seg.y;
% 
%                 xf1 = (R - delL) * cos(ths) + x;
%                 yf1 = (R - delL) * sin(ths) + y;
%             
%                 xf2 = (R + delL) * cos(ths) + x;
%                 yf2 = (R + delL) * sin(ths) + y;
%                 
%                 left_lane = [seg.init_pointL seg.last_pointL];
%                 right_lane = [seg.init_pointR seg.last_pointR];
% 
%                 plot(xf1,yf1,'k--');
%                 plot(xf2,yf2,'k--');
%                 plot(left_lane(1,:),left_lane(2,:),'bp')
%                 plot(right_lane(1,:),right_lane(2,:),'bp')
%             end
        end
        

    end
end

%% ============= Functions ============= 
function opt = FindBestIntvs(LP_l,LP_r,W_l,W_r,bnds,strs,order,line_segs,mode)
    
    % For initial/last segments, line_segs include only one line segment
    % information
    % For remaining segments, line_segs include two line segment
    % information
    [m,n] = size(strs);
    
    lb = (bnds(1)+1) * ones(1,n-1); ub = (bnds(end)-1) * ones(1,n-1);
    intcon = 1:(n-1); 
    delta = bnds(2)-bnds(1)-1;
    
    x0 = sort(bnds(1) + randperm(delta,n-1));
    
    A = zeros(n-2,n-1); b = -ones(n-2,1);
    for i=1:n-2
        A(i,i) = 1; A(i,i+1) = -1;
    end
%     options = optimoptions('surrogateopt','Display','off','UseParallel',true,'MaxFunctionEvaluations',500,'InitialPoints',x0);
    options = optimoptions('surrogateopt','Display','off','UseParallel',true,'MaxFunctionEvaluations',500);
    opt = struct();
    opt.plot_flag = false;
    opt.bnds = surrogateopt(@cost_func,lb,ub,intcon,A,b,[],[],options);
    
%     opt.plot_flag = false;
    cost_func(opt.bnds);
%     cost_func(x0);
%     disp(['Cost Func Val: ',num2str()])
    disp('-Error status for every segment- ')
    
    full_bnds = sort([bnds(1)-1 bnds(end) opt.bnds]);
%     full_bnds = sort([bnds(1)-1 bnds(end) x0]);
    for i=1:length(opt.err_maxs)
        disp(['Segment ',num2str(i),...
              ', Idx ',num2str(full_bnds(i)+1),...
              '~',num2str(full_bnds(i+1)),...
              '  Error: ',num2str(opt.err_maxs(i)),'m'])
        disp(strcat('Segment Type: ',strs(opt.f_min_idx,i)))
    end

    function f = cost_func(x) 
        fs = zeros(m,n);
        errs = zeros(m,n);
        
        for j=1:m
            % Each String Sequence
%             figure(j+1); hold on; grid on; axis equal;
            
            ord = order(j,:);
            
            segs = {};
            for k=1:n
                segs = [segs {0}];
                segs{k} = {};
            end
            
            switch mode
                case 'init'
                    segs = [segs line_segs];
                case 'mid'
                    segs = [segs line_segs(2) line_segs(1)];
                case 'last'
                    segs = [line_segs segs];
            end
            
            for k=1:n
                % Inside Each String Sequence
                seq_idx = find(ord == k);
                if seq_idx == 1
                    lb_ = bnds(1);
                    ub_ = x(seq_idx);
                    
                elseif seq_idx == n
                    lb_ = x(seq_idx-1)+1;
                    ub_ = bnds(end);
                else
                    lb_ = x(seq_idx-1)+1;
                    ub_ = x(seq_idx);
                end

                X_l = LP_l(1,lb_:ub_); X_r = LP_r(1,lb_:ub_);
                Y_l = LP_l(2,lb_:ub_); Y_r = LP_r(2,lb_:ub_);
                w_l = W_l(lb_:ub_); w_r = W_r(lb_:ub_);
                
                switch strs(j,seq_idx)
                    case "Front"
                        if isempty(getArrVal(segs,seq_idx-1))
                            error('Wrong Ordering Detected!')
                        else
                            seg1 = getArrVal(segs,seq_idx-1);
                            fixed_points = [seg1.last_pointL seg1.last_pointR];
                            [res,err] = CircleFitV3(X_l,X_r,Y_l,Y_r,w_l,w_r,...
                                                    opt.plot_flag,strs(j,seq_idx),...
                                                    fixed_points,true); 
                             
                        end

                    case "Back"
                        if isempty(getArrVal(segs,seq_idx+1))
                            error('Wrong Ordering Detected!')
                        else
                            seg2 = getArrVal(segs,seq_idx+1);
                            fixed_points = [seg2.init_pointL seg2.init_pointR];
                            [res,err] = CircleFitV3(X_l,X_r,Y_l,Y_r,w_l,w_r,...
                                                    opt.plot_flag,strs(j,seq_idx),...
                                                    fixed_points,true); 
%                             disp(['Back: ',num2str(err.tot)])
                        end
                        
                    case "Both"
                        if isempty(getArrVal(segs,seq_idx-1)) || isempty(getArrVal(segs,seq_idx+1))
                            error('Wrong Ordering Detected!')
                        else
                            seg1 = getArrVal(segs,seq_idx-1);
                            seg2 = getArrVal(segs,seq_idx+1);          
                            fixed_points = [seg1.last_pointL seg1.last_pointR seg2.init_pointL seg2.init_pointR];
                                                        
                            [res,err] = CircleFitV3(X_l,X_r,Y_l,Y_r,w_l,w_r,...
                                                    opt.plot_flag,strs(j,seq_idx),...
                                                    fixed_points,true); 
                            
%                             disp(['Both: ',num2str(err.tot)])
                        end

                    case "Free"
                        [res,err] = CircleFitV3(X_l,X_r,Y_l,Y_r,w_l,w_r,...
                                                opt.plot_flag,strs(seq_idx),...
                                                [],true); 
%                         disp(['Free: ',num2str(err.tot)])
                end

                segs{seq_idx} = res;
                fs(j,k) = err.tot; 
                errs(j,k) = err.max;
                
            end
            hold off;
        end
%         disp(fs)
        [f_,opt.f_min_idx] = min(sum(fs,2));
        opt.err_maxs = errs(opt.f_min_idx,:);
        f = f_;
    end
end

function val = getArrVal(arr,idx)
    if idx == 0
        val = arr{end};
    else
        val = arr{idx};
    end
end
