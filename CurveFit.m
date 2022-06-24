classdef CurveFit < handle
    properties
        plot_num = 1 % Plot Number
        Optimizer
        opt = struct()
        
    end
    methods
        %% Constructor
        function obj = CurveFit(Optimizer)
            obj.Optimizer = Optimizer;
        end
        
        %% Part 1: Curvature Segmentation using Dynamic Programming
        function obj = curvseg(obj, plot_flag)
            LP_l = reshape(obj.Optimizer.opt.lml(:,1),2,[]);
            n = size(LP_l,2);
            curv_l = zeros(obj.Optimizer.prev_num,n-2);
            curv_r = curv_l;
            
            disp('-Phase 1: Computing Discrete Curvature')
            for j=1:prev_num
                LP_l = reshape(obj.Optimizer.opt.lml(:,j),2,[]);
                LP_r = reshape(obj.Optimizer.opt.lmr(:,j),2,[]);
                n = size(LP_l,2);
            
                for i=2:n-1
                
                    p_im1_l = LP_l(:,i-1);
                    p_i_l = LP_l(:,i);
                    p_ip1_l = LP_l(:,i+1);
                
                    p_im1_r = LP_r(:,i-1);
                    p_i_r = LP_r(:,i);
                    p_ip1_r = LP_r(:,i+1);
                    
                    D_l = horzcat(p_i_l-p_im1_l,p_ip1_l-p_i_l);
                    D_r = horzcat(p_i_r-p_im1_r,p_ip1_r-p_i_r);
                    curv_l(j,i-1) = 2 * det(D_l) / (norm(p_i_l - p_im1_l) * norm(p_ip1_l - p_i_l) * norm(p_ip1_l - p_im1_l));
                    curv_r(j,i-1) = 2 * det(D_r) / (norm(p_i_r - p_im1_r) * norm(p_ip1_r - p_i_r) * norm(p_ip1_r - p_im1_r));
                    
                end
            end
            
            for j=1:prev_num
                curv_l(j,:) = wdenoise(curv_l(j,:),2,DenoisingMethod="Minimax");
                curv_r(j,:) = wdenoise(curv_r(j,:),2,DenoisingMethod="Minimax");
            end
            
            disp('-Phase 2: Performing Dynamic Programming for Curvature based Segmentation')
            [~, ~, obj.opt.intvs] = DP(curv_l(1,:),curv_r(1,:));
            % May need to run this via other programming language using GPU
            % --> To be done in the future

            n_intvs = length(obj.opt.intvs);
            obj.opt.kappas = zeros(1,n_intvs-1);

            for i=1:n_intvs-1
                lb = intvs_init(i); ub = intvs_init(i+1);
                kappa_mean_l = mean(curv_l(1,lb:ub));
                kappa_mean_r = mean(curv_r(1,lb:ub));
                obj.opt.kappas(i) = mean([kappa_mean_l; kappa_mean_r]);
            end

            if plot_flag
                figure(obj.plot_num);
                obj.plot_num = obj.plot_num + 1;

                subplot(2,1,1)
                plot(curv_l(1,:)); hold on; grid on;
                
                for i=1:n_intvs-1
                    lb = obj.opt.intvs(i); ub = obj.opt.intvs(i+1);
                    kappa_mean_l = mean(curv_l(1,lb:ub));
                    
                    plot([lb, ub],[kappa_mean_l, kappa_mean_l],'k--');
                    plot([lb, ub],[kappa_mean_l, kappa_mean_l],'bo');
                end
                
                ylabel('Curvature');
                title('Curvature based segmentation of 0m Previewed Left Lane')
                xlim([0,n]); ylim([-1.5e-3 1.5e-3])
                
                subplot(2,1,2);
                
                plot(curv_r(1,:)); hold on; grid on;
                
                for i=1:n_intvs-1
                    lb = obj.opt.intvs(i); ub = obj.opt.intvs(i+1);
                    kappa_mean_r = mean(curv_r(1,lb:ub));
                    
                    plot([lb, ub],[kappa_mean_r, kappa_mean_r],'k--');
                    plot([lb, ub],[kappa_mean_r, kappa_mean_r],'bo');
                end
                xlabel('Index'); ylabel('Curvature')
                title('Curvature based segmentation of 0m Previewed Right Lane')
                xlim([0,n]); ylim([-1.5e-3 1.5e-3])

            end

        end

        %% Part2: Optimize
        
        %% Phase 1: Find and parametrize straight lines
        function obj = optimizePh1(obj)
            obj.opt.line_segments = {};
            obj.opt.arc_segments = {};

             
            disp('-Phase 1: Identifying and Parametrization of Straight Lane Segments')
            thres = 1e-4;
%             acc_p = 0.1;
            n = length(obj.opt.kappas);
            obj.opt.st_intvs = [];

            for i=1:n
                if abs(obj.opt.kappas(i)) < thres
                    
                    lb_intv = 10*(obj.opt.intvs(i)-1)+1;
                    ub_intv = 10*(obj.opt.intvs(i+1)-1)+1;
                    lb_idx = find(obj.Optimizer.opt.reordered_lml_pc(3,:) == lb_intv);
                    ub_idx = find(obj.Optimizer.opt.reordered_lmr_pc(3,:) == ub_intv);

                    obj.opt.st_intvs = [obj.opt.st_intvs; lb_idx ub_idx];
                end
            end

            n = size(obj.opt.st_intvs,1);
            
            LP_l = obj.Optimizer.opt.reordered_lml_pc(1:2,:);
            LP_r = obj.Optimizer.opt.reordered_lmr_pc(1:2,:);
            
            w_l = obj.Optimizer.opt.w_l;
            w_r = obj.Optimizer.opt.w_r;
            
%             figure(2);
%             obj.plot_num = obj.plot_num+1;
%             plot(LP_l(1,:),LP_l(2,:),'r.'); hold on; grid on; axis equal;
%             plot(LP_r(1,:),LP_r(2,:),'g.');
%             

            for i=1:n
                lb_idx = obj.opt.st_intvs(i,1); 
                ub_idx = obj.opt.st_intvs(i,2);
                x0 = [lb_idx, ub_idx];
                [lb, ub] = findBestIntv(LP_l,LP_r,w_l,w_r,x0);

                X_l = LP_l(1,lb:ub); X_r = LP_r(1,lb:ub);
                Y_l = LP_l(2,lb:ub); Y_r = LP_r(2,lb:ub);

                disp(['Line Segment ID: ', num2str(i),', Segment Index: ',num2str(lb),' ~ ',num2str(ub)])
                [res, err] = LineFitV2(X_l,X_r,Y_l,Y_r,...
                                       w_l(lb:ub),w_r(lb:ub),false);
                seg = struct();
                seg.res = res;
                seg.err = err;
                seg.status = 'line';
                seg.bnds = [lb, ub];
                obj.opt.line_segments = [obj.opt.line_segments {seg}];
                obj.opt.st_intvs(i,:) = seg.bnds;
            end
        end

        %% Phase 2: Fill in remaining parts with circular arcs
        function obj = optimizePh2(obj)
            
            LP_l = obj.Optimizer.opt.reordered_lml_pc(1:2,:);
            LP_r = obj.Optimizer.opt.reordered_lmr_pc(1:2,:);
            
            w_l = obj.Optimizer.opt.w_l;
            w_r = obj.Optimizer.opt.w_r;

            disp('-Phase 2: Filling in Remaining Parts with Circular Arcs')
            
            n = size(obj.opt.st_intvs,1);
            obj.opt.arc_intvs = zeros(n+1,2);
            obj.opt.arc_intvs(1,1) = 1;
            obj.opt.arc_intvs(1,2) = obj.opt.st_intvs(1,1)-1;

            for i=1:n-1
                obj.opt.arc_intvs(i+1,1) = obj.opt.st_intvs(i,2)+1;
                obj.opt.arc_intvs(i+1,2) = obj.opt.st_intvs(i+1,1)-1;
            end

            obj.opt.arc_intvs(end,1) = obj.opt.st_intvs(end,2)+1;
            obj.opt.arc_intvs(end,2) = size(obj.Optimizer.opt.reordered_lml_pc,2);
            
            %% One-way optimization for first and last segments
            disp('-Arc Approximation for first segment-')
            
            lb = obj.opt.arc_intvs(1,2)-300;
            ub = obj.opt.arc_intvs(1,2)-1;
            seg = obj.opt.line_segments{1};
            
            fixed_status = 'Back';
            adjacent_seg_type = 'line';
            fixed_points = [seg.res.init_pointL seg.res.init_pointR];
            prev_D = 0;
            th = seg.res.theta;

            while true
%                 disp(['Current Lower Bound: ',num2str(lb)])
                X_l = LP_l(1,lb:ub); Y_l = LP_l(2,lb:ub);
                X_r = LP_r(1,lb:ub); Y_r = LP_r(2,lb:ub);
                W_l = w_l(lb:ub); W_r = w_r(lb:ub);

                [res, err] = CircleFitV3(X_l,X_r,Y_l,Y_r,W_l,W_r,...
                                         false,fixed_status,adjacent_seg_type,...
                                         fixed_points,prev_D,th);

                err_tot = [err.full_l err.full_r];
                
                if length(find(err_tot >= 10*1e-2)) >= 3
                    disp('-Threshold for arc approximation exceeded')
                    disp(['-Segment: Idx ',num2str(lb+1),' ~ ',num2str(ub)])

                    seg = struct();
                    seg.res = prev_res;
                    seg.err = prev_err;
                    seg.status = 'init'; % Marking initial segment approximation
                    seg.bnds = [lb+1 ub];

                    obj.opt.arc_segments = [{seg} obj.opt.arc_segments];

                    adjacent_seg_type = 'arc';
                    fixed_points = [prev_res.init_pointL prev_res.init_pointR];
                    prev_D = prev_res.D;
                    th = prev_res.th_init;
                    
                    ub = lb;
                    lb = lb - 300;

                    if lb <= 0
                        disp('Arc Approximation for first segment finished')
                        lb = 1;
                        X_l = LP_l(1,lb:ub); Y_l = LP_l(2,lb:ub);
                        X_r = LP_r(1,lb:ub); Y_r = LP_r(2,lb:ub);
                        W_l = w_l(lb:ub); W_r = w_r(lb:ub);

                        [res, err] = CircleFitV3(X_l,X_r,Y_l,Y_r,W_l,W_r,...
                                               false,fixed_status,adjacent_seg_type,...
                                               fixed_points,prev_D,th);
                        seg = struct();
                        seg.res = res;
                        seg.err = err;
                        seg.status = 'init';
                        seg.bnds = [1 ub];
                        
                        obj.opt.arc_segments = [{seg} obj.opt.arc_segments];
                        break;
                    end
                else
                    prev_res = res;
                    prev_err = err;
                    lb = lb-1;
                    
                    if lb == 0
                        disp('Arc Approximation for first segment finished')
                        obj.opt.arc_segments = [{prev_res} obj.opt.arc_segments];
                        break;
                    end
                end
            end

            disp('-Arc Approximation for last segment-')
            lb = obj.opt.arc_intvs(end,2)+1;
            ub = obj.opt.arc_intvs(end,2)+300;
            n = size(obj.Optimizer.reordered_lml_pc,2);
            seg = obj.opt.line_segments{end};
            
            fixed_status = 'Front';
            adjacent_seg_type = 'line';
            fixed_points = [seg.res.last_pointL seg.res.last_pointR];
            prev_D = 0;
            th = seg.res.theta;
        
            while true
%                 disp(['Current Lower Bound: ',num2str(lb)])
                X_l = LP_l(1,lb:ub); Y_l = LP_l(2,lb:ub);
                X_r = LP_r(1,lb:ub); Y_r = LP_r(2,lb:ub);
                W_l = w_l(lb:ub); W_r = w_r(lb:ub);

                [res, err] = CircleFitV3(X_l,X_r,Y_l,Y_r,W_l,W_r,...
                                         false,fixed_status,adjacent_seg_type,...
                                         fixed_points,prev_D,th);

                err_tot = [err.full_l err.full_r];
                
                if length(find(err_tot >= 10*1e-2)) >= 3
                    disp('-Threshold for arc approximation exceeded')
                    disp(['-Segment: Idx ',num2str(lb),' ~ ',num2str(ub-1)])

                    seg = struct();
                    seg.res = prev_res;
                    seg.err = prev_err;
                    seg.status = 'last'; % Marking initial segment approximation
                    seg.bnds = [lb ub-1];

                    obj.opt.arc_segments = [{seg} obj.opt.arc_segments];

                    adjacent_seg_type = 'arc';
                    fixed_points = [prev_res.init_pointL prev_res.init_pointR];
                    prev_D = prev_res.D;
                    th = prev_res.th_init;
                    
                    lb = ub;
                    ub = lb + 300;

                    if ub > n
                        disp('Arc Approximation for last segment finished')
                        lb = 1;
                        X_l = LP_l(1,lb:ub); Y_l = LP_l(2,lb:ub);
                        X_r = LP_r(1,lb:ub); Y_r = LP_r(2,lb:ub);
                        W_l = w_l(lb:ub); W_r = w_r(lb:ub);

                        [res, err] = CircleFitV3(X_l,X_r,Y_l,Y_r,W_l,W_r,...
                                               false,fixed_status,adjacent_seg_type,...
                                               fixed_points,prev_D,th);
                        seg = struct();
                        seg.res = res;
                        seg.err = err;
                        seg.status = 'last';
                        seg.bnds = [lb n];
                        
                        obj.opt.arc_segments = [{seg} obj.opt.arc_segments];
                        break;
                    end
                else
                    prev_res = res;
                    prev_err = err;
                    ub = ub+1;
                    
                    if ub > n
                        disp('Arc Approximation for first segment finished')
                        obj.opt.arc_segments = [{prev_res} obj.opt.arc_segments];
                        break;
                    end
                end
            end


            %% Two-way optimization for segments in between
            disp('-Arc Approximation for remaining segment-')

            
        end

        %% Phase 3: Shift Segments to match lateral position (Optimization)
        function obj = optimizePh3(obj)
        end

        function plotRes(obj)
            line_segs = obj.opt.line_segments;
            arc_segs = obj.opt.arc_segments;
            
            LP_l = obj.Optimizer.opt.reordered_lml_pc(1:2,:);
            LP_r = obj.Optimizer.opt.reordered_lmr_pc(1:2,:);

            % Plot Line Segments 
            figure(2); hold on; axis equal; grid on;
            plot(LP_l(1,:),LP_l(2,:),'r.');
            plot(LP_r(1,:),LP_r(2,:),'g.');

            n = length(line_segs);
            for i=1:n
                seg = line_segs{i};
                left_lane = [seg.res.init_pointL seg.res.last_pointL];
                right_lane = [seg.res.init_pointR seg.res.last_pointR];
              
                plot(left_lane(1,:),left_lane(2,:),'k--')
                plot(right_lane(1,:),right_lane(2,:),'k--')
                plot(left_lane(1,:),left_lane(2,:),'bp')
                plot(right_lane(1,:),right_lane(2,:),'bp')
            end

            % Plot Arc Segments
            n = length(arc_segs);
            for i=1:n
                seg = arc_segs{i};
                ths = linspace(seg.res.th_init, seg.res.th_last, 1e4);
                
                R = seg.res.R; delL = seg.res.delL;
                x = seg.res.x; y = seg.res.y;

                xf1 = (R - delL) * cos(ths) + x;
                yf1 = (R - delL) * sin(ths) + y;
            
                xf2 = (R + delL) * cos(ths) + x;
                yf2 = (R + delL) * sin(ths) + y;
                
                left_lane = [seg.res.init_pointL seg.res.last_pointL];
                right_lane = [seg.res.init_pointR seg.res.last_pointR];

                plot(xf1,yf1,'k--');
                plot(xf2,yf2,'k--');
                plot(left_lane(1,:),left_lane(2,:),'bp')
                plot(right_lane(1,:),right_lane(2,:),'bp')
            end
        end
        

    end
end

%% ============= Functions ============= 
%% Dynamic Programming for Curvature Segmentation
function [M, sep, intvs] = DP(curv_l, curv_r)
    n = length(curv_l);
    M = zeros(n,n);
    sep = zeros(n,n);
    E_cost = 0.05;

    for i=1:n-1
        M(i,i+1) = StepFit(i,i+1) + E_cost;
    end

    for i=3:n
        disp(['Iteration: ',num2str(i)])
        for j=1:n-i+1
            lst = StepFit(j,i+j-1) + E_cost;

            for k=j+1:i+j-2
                lst = [lst M(j,k) + M(k,i+j-1)];
            end
            
            [M(j,i+j-1), idx] = min(lst);

            if idx == 1
                sep(j,i+j-1) = 0;
            else
                sep(j,i+j-1) = j + idx-1;
            end
            
        end
    end
    
    intvs = sort([1,n,backTrack(sep,1,n)]);
    
    function intvs = backTrack(sep, search_row, search_col)
        intvs = [];
        if sep(search_row, search_col) ~= 0
            search_idx = sep(search_row, search_col);
            intvs = [intvs,...
                     search_idx, ...
                     backTrack(sep,search_row,search_idx),...
                     backTrack(sep,search_idx,search_col)];
        end

    end


    function err = StepFit(idx1, idx2)
        kappa = mean([curv_l(idx1:idx2) curv_r(idx1:idx2)]);
        err = sum(abs([curv_l(idx1:idx2) curv_r(idx1:idx2)] - kappa));
    end
end

%% Find Best Interval for Line Fitting using Surrogate Optimization
function [lb, ub] = findBestIntv(LP_l,LP_r,w_l,w_r,x0)
    options = optimoptions('surrogateopt','Display','none','UseParallel',true,'MaxFunctionEvaluations',200);

    opt = surrogateopt(@costfunc,[x0(1) x0(1)],[x0(2) x0(2)],[1 2],[1 -1],-100,[],[],options);
    lb = opt(1); ub = opt(2);

    %% Cost Function for Surrogate Optimization
    function f = costfunc(x)
        lb_ = x(1); ub_ = x(2);
        X_l = LP_l(1,lb_:ub_);
        Y_l = LP_l(2,lb_:ub_);
        X_r = LP_r(1,lb_:ub_);
        Y_r = LP_r(2,lb_:ub_);
        [~,err] = LineFitV2(X_l,X_r,Y_l,Y_r,...
                              w_l(lb_:ub_), w_r(lb_:ub_),false);

        f.Fval = 1/(ub_ - lb_ + 1)^(3.4) * err.tot;
    end
    
end


