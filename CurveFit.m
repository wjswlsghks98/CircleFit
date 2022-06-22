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
        function obj = optimize(obj)
            obj.opt.line_segments = {};
            obj.opt.arc_segments = {};

            % Phase 1: Find and parametrize straight lines 
            disp('-Phase 1: Identifying and Parametrization of Straight Lane Segments')
            thres = 1e-4;
%             acc_p = 0.1;
            acc_thres = 10*1e-2;
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
            
%             w_l = obj.Optimizer.opt.w_l;
%             w_r = obj.Optimizer.opt.w_r;
            w_l = ones(1,size(LP_l,2));
            w_r = w_l;
            
            figure(2);
            obj.plot_num = obj.plot_num+1;
            plot(LP_l(1,:),LP_l(2,:),'r.'); hold on; grid on; axis equal;
            plot(LP_r(1,:),LP_r(2,:),'g.');
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
                                                       w_l(lb:ub),w_r(lb:ub),true);
                seg = struct();
                seg.res = res;
                seg.err = err;
                seg.bnds = [lb, ub];
                obj.opt.line_segments = [obj.opt.line_segments {seg}];
            end
            
            % Phase 2: Fill in remaining parts with circular arcs
            disp('-Phase 2: Filling in Remaining Parts with Circular Arcs')
            
            
            
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
        [res,err] = LineFitV2(X_l,X_r,Y_l,Y_r,...
                              w_l(lb_:ub_), w_r(lb_:ub_),false);

        f.Fval = 1/(res.L)^(3.6) * err.tot;
    end
    
end

