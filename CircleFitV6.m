classdef CircleFitV6 < handle
% [CircleFitV6: Nonlinear Least Squares based Multi-arc spline]
% Using curvature and arc_length as parameters
% Curve is continuously connected (2D and heading angle continuity is guaranteed)
% [Optimization Step]
% --Initial Step --
% 1. Given lane points, perform recursive line interpolation to find
% polyline segment intervals
% 2. Merge line intervals for appropriate circular fitting
% -- Data Association --
% 3. Using each segment information, classify and divide lane points.
% -- NLS optimization
% 4. Perform NLS optimization with fixed data association obtained at step 
% 3 and check if all segments are valid.
% * If all segments are valid, end optimization
% * If not all segments are valid, go back to step 3, and add iteration
% counter.
% * If iteration number exceeded maximum iteration limit but the
% parametrization is still invalid, add segment(s) at invalid segment(s)

% Implemented by Jin Hwan Jeon, 2022
% email: jordan98@kaist.ac.kr
    properties
        LP % Data points
        cov % Data point covariance 
        w % Data point 1D weight, lateral std^2
        thres % Chi-Square test accuracy threshold: 0.9, 0.95, 0.99
        intv % Data point interval of interest
        num_seg % Number of target arc segments
        init_segments = {} % Initial arc segments
        opt = struct() % Optimized Results
        segments = {} % Each arc segment parameters        
        fig_num % Figure number
    end

    methods
        %% Constructor
        function obj = CircleFitV6(LP,w,cov,thres,intv,fig_num)
            obj.LP = LP(1:2,:);
            obj.w = w;
            obj.cov = cov;
            obj.thres = thres;
            obj.intv = intv;     
            obj.fig_num = fig_num;
        end

        %% Initialize
        function obj = initialize(obj)           
            % [Part 1: Create initial value segments]
            % Perform divide-conquer based recursive line segmentation for
            % finding initial index intervals. Then merge intervals that
            % can be stably connected by circluar fitting.
            % Initial number of segments is determined by this process
            
            obj.findInitialSegments();                     

            % [Part 2: Create initial value for optimization]
            obj.opt.X0 = zeros(3+2*obj.num_seg,1);
            obj.opt.X0(1:2) = obj.LP(:,obj.intv(1));
            obj.opt.X0(3) = atan2(obj.LP(2,obj.intv(1)+1) - obj.LP(2,obj.intv(1)),...
                                  obj.LP(1,obj.intv(1)+1) - obj.LP(1,obj.intv(1)));
            for i=1:obj.num_seg
                obj.opt.X0(3+2*i-1:3+2*i) = [obj.init_segments{i}.kappa;
                                             obj.init_segments{i}.L];
            end
            obj.opt.init_num_seg = obj.num_seg;
            obj.opt.initX0 = obj.opt.X0;
            
            % [Part 3: Data association]
            modifiedLP = [obj.LP; obj.intv(1):obj.intv(end)];
            obj.opt.assoc = classifyLP(modifiedLP,obj.init_segments,true);
        end
        
        %% Find Initial Segments
        function obj = findInitialSegments(obj)
            % There are total of two steps in this function.
            % [Step 1: Line based recursive segmentation]            

            disp('=========[Initial Segmentation]=========');
            n = size(obj.LP,2);
            obj.opt.initL_intvs = unique(obj.segmentation([1,n],1));

            % [Step 2: Merge line intervals for valid circular approximation]

            obj.mergeSeg(obj.opt.initL_intvs);
            disp('=========[Initial Segmentation Finished]=========');
        end

        %% Optimize
        function obj = optimize(obj)
            options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',false,...
                                   'Display','off','Algorithm','trust-region-reflective',...
                                   'UseParallel',true,'MaxFunctionEvaluations',inf,...
                                   'StepTolerance',1e-6,'MaxIterations',inf,'CheckGradients',false);
            max_iter = 10; 
            obj.opt.valid = false;
            
            while ~obj.opt.valid
                iter = 1;
                disp(['[Batch Optimization Starts, number of segments: ',num2str(obj.num_seg),']'])
                while iter <= max_iter
                    % Need to set lower bounds for arc lengths
                    lb1 = [-inf,-inf,-inf];
                    lb2 = [-inf 20];
                    lb = horzcat(lb1,repmat(lb2,1,obj.num_seg));
                    
                    tic;
                    [obj.opt.fullx,~,obj.opt.res,~,~,~,obj.opt.jac] = lsqnonlin(@obj.cost_func,obj.opt.X0,lb,[],options);
                    obj.opt.time = toc;
%                     disp(['Time taken for NLS optimization: ',num2str(obj.opt.time),'s'])
                    % Post process optimization results to check validity of
                    % arc spline approximation
                    obj.opt.segments = obj.segments;
                    m = length(obj.segments);
                    x0 = obj.opt.fullx(1); y0 = obj.opt.fullx(2); tau = obj.opt.fullx(3);
                    kappa = obj.opt.fullx(4); L = obj.opt.fullx(5);
        
                    xc = x0 - 1/kappa * sin(tau); yc = y0 + 1/kappa * cos(tau);
                    init_point = [x0; y0];
                    last_point = [xc + 1/kappa * sin(tau+kappa*L);
                                  yc - 1/kappa * cos(tau+kappa*L)];
                    obj.opt.segments{1}.xc = xc;
                    obj.opt.segments{1}.yc = yc;
                    obj.opt.segments{1}.kappa = kappa;
                    obj.opt.segments{1}.L = L;
                    obj.opt.segments{1}.init_point = init_point;
                    obj.opt.segments{1}.last_point = last_point;
                    obj.opt.segments{1}.th_init = atan2(init_point(2)-yc,init_point(1)-xc);
                    obj.opt.segments{1}.th_last = atan2(last_point(2)-yc,last_point(1)-xc);            
                    
                    for i=2:m                                   
                        tau = tau + kappa * L;
                        xc = xc + 1/kappa * sin(tau);
                        yc = yc - 1/kappa * cos(tau);
                        kappa = obj.opt.fullx(3+2*i-1);
        
                        xc = xc - 1/kappa * sin(tau);
                        yc = yc + 1/kappa * cos(tau);
                        L = obj.opt.fullx(3+2*i);
        
                        init_point = last_point;
                        last_point = [xc + 1/kappa * sin(tau+kappa*L);
                                      yc - 1/kappa * cos(tau+kappa*L)];
                        
                        obj.opt.segments{i}.xc = xc;
                        obj.opt.segments{i}.yc = yc;
                        obj.opt.segments{i}.kappa = kappa;
                        obj.opt.segments{i}.L = L;
                        obj.opt.segments{i}.init_point = init_point;
                        obj.opt.segments{i}.last_point = last_point;
                        obj.opt.segments{i}.th_init = atan2(init_point(2)-yc,init_point(1)-xc);
                        obj.opt.segments{i}.th_last = atan2(last_point(2)-yc,last_point(1)-xc);
                    end
                    
                    % Check if all arc segments are valid
                    obj.computeError(false);    
                    
                    obj.plotRes();

                    % If all 
                    if obj.opt.valid      
                        obj.computeError(true);
                        break;
                    else
                        modifiedLP = [obj.LP; obj.intv(1):obj.intv(2)];
                        obj.opt.assoc = classifyLP(modifiedLP,obj.opt.segments,false);
                        obj.opt.X0 = obj.opt.fullx;
                        iter = iter + 1;
                    end
                end
                
                

                if obj.opt.valid
                    break;
                else
                    % Increase segment
                    obj.createNewX0();
                end
            end
        end

        %% Cost Function
        function res = cost_func(obj,x0)
            obj.segments = {};
            % [Part 1: Precompute parameter jacobians and calculate circle centers]
            obj.precompJac(x0);  
            
            % [Part 2: Compute full jacobian and residual]
            % Part 2-1: Normal Measurement Jacobian using chain rule
            obj.CreateMEBlock();
            % Part 2-2: Initial point and Last point Measurement Jacobian
            obj.CreateMEBlock2(x0);

            obj.Merge();
            res = obj.opt.res;                                            
        end
        
        %% Precompute Jacobians
        function obj = precompJac(obj,x0)
            tau = x0(3); kappa = x0(4); L = x0(5);
            xc = x0(1) - 1/kappa * sin(tau); yc = x0(2) + 1/kappa * cos(tau);
            init_point = x0(1:2);
            last_point = init_point + 1/kappa * [sin(tau+kappa*L)-sin(tau); -cos(tau+kappa*L)+cos(tau)];
            jac_ = zeros(3,3+2*obj.num_seg);
            jac_(:,1:4) = [1, 0, -1/kappa * cos(tau), 1/kappa^2 * sin(tau);
                           0, 1, -1/kappa * sin(tau), -1/kappa^2 * cos(tau);
                           0, 0,                   0,                     1];
            seg = struct();
            seg.xc = xc; seg.yc = yc; seg.kappa = kappa;
            seg.Cjac = jac_; % Chain Rule Based Pre-computed Jacobian
            seg.init_point = init_point;
            seg.last_point = last_point;
            seg.th_init = atan2(init_point(2) - yc,init_point(1) - xc);
            seg.th_last = atan2(last_point(2) - yc,last_point(1) - xc);
            obj.segments = [obj.segments {seg}];

            for i=2:obj.num_seg                               
                % Propagate Variables and Jacobians
                tau = tau + kappa * L;
                xc = xc + (1/kappa) * sin(tau);
                yc = yc - (1/kappa) * cos(tau);
                delta_jac = zeros(3,3+2*obj.num_seg);
                delta_jac(1,3) = 1/kappa * cos(tau);
                delta_jac(2,3) = 1/kappa * sin(tau);
                
                for j=1:i
                    Lj = x0(3+2*j); kappaj =x0(3+2*j-1);

                    if j == i-1
                        delta_jac(1,3+2*j-1) = L/kappa * cos(tau) - 1/kappa^2 * sin(tau);
                        delta_jac(2,3+2*j-1) = 1/kappa^2 * cos(tau) + L/kappa * sin(tau);
                        delta_jac(3,3+2*j-1) = -1; 
                    elseif j == i
                        delta_jac(3,3+2*j-1) = 1; 
                    else                        
                        delta_jac(1,3+2*j-1) = Lj/kappa * cos(tau);
                        delta_jac(2,3+2*j-1) = Lj/kappa * sin(tau);                        
                    end

                    if j~=i
                        delta_jac(1,3+2*j) = kappaj/kappa * cos(tau);
                        delta_jac(2,3+2*j) = kappaj/kappa * sin(tau);
                    end
                end

                kappa = x0(3+2*i-1); % Ki --> Ki+1
                
                xc = xc - (1/kappa) * sin(tau);
                yc = yc + (1/kappa) * cos(tau);
                % Xc_(i+1) = Xc_(i) + (1/Ki - 1/Ki+1) * sin(tau0 + K1L1 + ... KiLi);
                % Yc_(i+1) = Yc_(i) - (1/Ki - 1/Ki+1) * cos(tau0 + K1L1 + ... KiLi);

                delta_jac(1,3) = delta_jac(1,3) - 1/kappa * cos(tau);
                delta_jac(2,3) = delta_jac(2,3) - 1/kappa * sin(tau);
                for j=1:i
                    Lj = x0(3+2*j); kappaj =x0(3+2*j-1);    

                    if j == i-1
                        delta_jac(1,3+2*j-1) = delta_jac(1,3+2*j-1) - L/kappa * cos(tau);
                        delta_jac(2,3+2*j-1) = delta_jac(2,3+2*j-1) - L/kappa * sin(tau);
                    elseif j == i
                        delta_jac(1,3+2*j-1) = 1/kappa^2 * sin(tau);
                        delta_jac(2,3+2*j-1) = 1/kappa^2 * cos(tau);
                    else                        
                        delta_jac(1,3+2*j-1) = delta_jac(1,3+2*j-1) - Lj/kappa * cos(tau);
                        delta_jac(2,3+2*j-1) = delta_jac(2,3+2*j-1) - Lj/kappa * sin(tau);                        
                    end

                    if j~=i
                        delta_jac(1,3+2*j) = delta_jac(1,3+2*j) - kappaj/kappa * cos(tau);
                        delta_jac(2,3+2*j) = delta_jac(2,3+2*j) - kappaj/kappa * sin(tau);
                    end
                end

                jac_ = jac_ + delta_jac;

                L = x0(3+2*i); % Li --> Li+1

                seg = struct();
                seg.xc = xc; seg.yc = yc; seg.kappa = kappa;
                seg.Cjac = jac_; % Chain Rule Based Pre-computed Jacobian
                init_point = last_point;
                last_point = init_point + 1/kappa * [sin(tau+kappa*L)-sin(tau); -cos(tau+kappa*L)+cos(tau)];
                seg.init_point = init_point;
                seg.last_point = last_point;
                seg.th_init = atan2(init_point(2) - yc,init_point(1) - xc);
                seg.th_last = atan2(last_point(2) - yc,last_point(1) - xc);
                obj.segments = [obj.segments {seg}];
            end
            
%             figure(1); hold on; grid on; axis equal;
%             m = length(obj.segments);
%             for i=1:m
%                 plot(obj.segments{i}.init_point(1),...
%                      obj.segments{i}.init_point(2),'kx');
%                 plot(obj.segments{i}.last_point(1),...
%                      obj.segments{i}.last_point(2),'rx');
%             end
%             plot(obj.LP(1,obj.intv(1):obj.intv(2)),...
%                  obj.LP(2,obj.intv(1):obj.intv(2)),'b.');

        end

        %% Create Normal Measurement Jacobian
        function obj = CreateMEBlock(obj)
            m = length(obj.segments);            
            
            obj.opt.ME_vec = zeros(obj.intv(end)-obj.intv(1)+1+(obj.num_seg-1),1);
            cnt = 1;
            for i=1:m
                idxs = obj.opt.assoc{i}.LP_idxs;
                n = length(idxs);
                xc = obj.segments{i}.xc;
                yc = obj.segments{i}.yc;
                kappa = obj.segments{i}.kappa;
                
                for j =1:n                    
                    cov_ = obj.w(idxs(j));
                    xk = obj.LP(1,idxs(j));
                    yk = obj.LP(2,idxs(j));
                    d = sqrt((xc - xk)^2 + (yc - yk)^2);
                    resi = d - abs(1/kappa);
                    obj.opt.ME_vec(cnt) = InvMahalanobis(resi,cov_);  
                    cnt = cnt + 1;
                end
            end                                 
        end

        %% Create Intial, Last Point Measurement Jacobian
        function obj = CreateMEBlock2(obj,x0)
            
            obj.opt.ME_vec2 = zeros(2*(obj.num_seg+1),1);
            % Initial Point
            cov1 = reshape(obj.cov(:,obj.intv(1)),2,2);
            obj.opt.ME_vec2(1:2) = InvMahalanobis(x0(1:2) - obj.LP(:,obj.intv(1)),cov1);
            
            % Final Point
            cov2 = reshape(obj.cov(:,obj.intv(end)),2,2);            
            kappas = x0(3+(1:2:2*obj.num_seg))';
            Ls = x0(3+(2:2:2*obj.num_seg))';
            alpha = x0(3) + kappas * Ls';

            % compute residual
            xc = obj.segments{end}.xc;
            yc = obj.segments{end}.yc;
            Z_pred = [xc + 1/kappas(end) * sin(alpha);
                      yc - 1/kappas(end) * cos(alpha)];
            
            obj.opt.ME_vec2(3:4) = InvMahalanobis(Z_pred - obj.LP(:,obj.intv(end)),cov2);
        end

        %% Merge
        function obj = Merge(obj)
            obj.opt.res = vertcat(obj.opt.ME_vec,...
                                  obj.opt.ME_vec2);
        end
        
        %% Compute Error for every segment
        function obj = computeError(obj,disp_flag)
            m = length(obj.opt.segments);
            obj.opt.valid = true;
            obj.opt.err = {};

%             disp('Error Analysis')
            for i=1:m                
                err_ = struct();
                idxs = obj.opt.assoc{i}.LP_idxs;
                n = length(idxs);
                err_.full = zeros(1,n);
                
                err_.wfull = zeros(1,n);
                err_.wthres = sqrt(chi2inv(obj.thres,2));
                err_.invalid_idxs = [];
                
                xc = obj.opt.segments{i}.xc;
                yc = obj.opt.segments{i}.yc;
                R = abs(1/obj.opt.segments{i}.kappa);
                
                cnt = 0;
                for j=1:n                   
                    lp = obj.LP(:,idxs(j));
                    th = atan2(lp(2)-yc,lp(1)-xc);
                    lp_pred = [xc + R * cos(th);
                               yc + R * sin(th)];
                    diff = lp_pred - lp;
                    cov_ = reshape(obj.cov(:,idxs(j)),2,2);
                    err_.full(j) = sqrt(diff' * diff);
                    err_.wfull(j) = sqrt(diff' / cov_ * diff);
                    
                    if err_.wfull(j) > err_.wthres
                        cnt = cnt + 1;
                        err_.invalid_idxs = [err_.invalid_idxs idxs(j)];
                    end
                end
                err_.emax = max(err_.full);
                
                if cnt < 3
                    if disp_flag
                        disp(['Segment ',num2str(i),': Valid, Maximum Error: ',num2str(err_.emax)])
                    end
                    err_.valid = true;
                else
                    
%                     disp(['Segment ',num2str(i),': Invalid, Maximum Error: ',num2str(err_.emax)])
                    err_.valid = false;
                    obj.opt.valid = false;
                end

                obj.opt.err = [obj.opt.err {err_}];
            end
        end

        %% Plot results
        function plotRes(obj)
            m = length(obj.segments);
            figure(obj.fig_num); 
            trig = false;
            for i=1:m
                xc = obj.opt.segments{i}.xc;
                yc = obj.opt.segments{i}.yc;                               
                
                R = abs(1/obj.opt.segments{i}.kappa);
                th_init = obj.opt.segments{i}.th_init;
                th_last = obj.opt.segments{i}.th_last;
                th = linspace(th_init,th_last,1e3);

                LPs = [xc + R * cos(th);
                       yc + R * sin(th)];
                init_point = obj.opt.segments{i}.init_point;
                last_point = obj.opt.segments{i}.last_point;
                point = [init_point last_point];
                
                p_approx = plot(LPs(1,:),LPs(2,:),'k--'); hold on; grid on; axis equal;

                idx = obj.opt.assoc{i}.LP_idxs;
                p_data = plot(obj.LP(1,idx),...
                              obj.LP(2,idx),'r.');
                p_cts = plot(point(1,:),point(2,:),'bp');

                if ~isempty(obj.opt.err{i}.invalid_idxs)
                    idxs = obj.opt.err{i}.invalid_idxs;
                    p_vio = plot(obj.LP(1,idxs),...
                                 obj.LP(2,idxs),'cs');
                    trig = true;
                end                
            end
            

            xlabel('Global X'); ylabel('Global Y'); 
            title(['Batch Arc Spline Fitting using ',num2str(obj.num_seg),' segments']);
            
            if trig
                legend([p_approx,p_data,p_cts,p_vio],...
                        'Piecewise Arc Splines','Data Points',...
                        'Control Points','Threshold Violated Points')
            else
                legend([p_approx,p_data,p_cts],...
                        'Piecewise Arc Splines','Data Points',...
                        'Control Points')
            end       
            hold off;
        end
        
        %% Create New X0 (initial value) for introducing new segment
        function obj = createNewX0(obj)
            n = length(obj.opt.segments);
            % Initial Point and initial heading angle is set equal
            
            x0 = [];
            cpyd_segments = {}; % To associate newly divided segments
            invalid_cnt = 0;
            for i=1:n
                if ~obj.opt.err{i}.valid % Current Segment is not valid
                    
                    cpyd_segments = [cpyd_segments obj.opt.segments{i} obj.opt.segments{i}];
                    cpyd_segments{end-1}.th_last = 1/2 * (obj.opt.segments{i}.th_last + obj.opt.segments{i}.th_init);
                    cpyd_segments{end}.th_init = 1/2 * (obj.opt.segments{i}.th_last + obj.opt.segments{i}.th_init);

                    % modify init_theta, last_theta ...
                    x0 = [x0; obj.opt.segments{i}.kappa; 1/2 * obj.opt.segments{i}.L;
                              obj.opt.segments{i}.kappa; 1/2 * obj.opt.segments{i}.L];    
                    invalid_cnt = invalid_cnt + 1;
                else
                    cpyd_segments = [cpyd_segments obj.opt.segments{i}];
                    x0 = [x0; obj.opt.segments{i}.kappa; obj.opt.segments{i}.L];
                end
            end
            modifiedLP = [obj.LP; obj.intv(1):obj.intv(end)];
            obj.opt.assoc = classifyLP(modifiedLP,cpyd_segments,false);
            obj.opt.X0 = [obj.opt.fullx(1:3); x0];  
            obj.num_seg = obj.num_seg + invalid_cnt;
        end
        
        %% Initial segmentation
        function intvs = segmentation(obj,intv,depth)
            % (1). Connect first and end intervals with a line
            % (2). Find the maximum vertical error data point
            % (3). Divide the data points of interest into two w.r.t point found at (b).
            % (4). Repeat (1)~(3) for every divided segments.
            % * If maximum error computed at (2) is lower than threshold, 
            % stop dividing current segment further and propagate backwards. 

            % Adjust this value to increase or decrease the number of line
            % interpolation. Typically, values smaller than 0.3m is
            % recommended.
            line_acc_thres = 0.2;

            %% Create Line and visualize current segmentation state 
            init_point = obj.LP(:,intv(1));
            last_point = obj.LP(:,intv(2));
            m = (last_point(2) - init_point(2))/(last_point(1) - init_point(1));
            n = -m*init_point(1) + init_point(2);
            % Visualize 
%             x = linspace(init_point(1),last_point(1),1e3);
%             y = m*x+n;
%             plot(x,y,'k--');
%             plot([x(1) x(end)],[y(1), y(end)],'bs')
            
            %% Find existence of intersection
            % Even if more than 1 intersection occurs, divide current segment with
            % two segments. For stability, search intersection from the middle
            X = obj.LP(1,intv(1):intv(2)); Y = obj.LP(2,intv(1):intv(2));
            d = (Y - m * X - n)/sqrt(m^2 + 1); % signed distance to line
            
            % If the maximum error caused by "Line Interpolation" is lower than
            % threshold, end segmentation(return current search interval). 
            % Else, divide segmentation problem using the maximum error point
            [max_val,max_idx] = max(abs(d));
%             plot(obj.LP(1,intv(1)+max_idx-1),obj.LP(2,intv(1)+max_idx-1),'cs');
%             pause(0.3);
            
            % End segmentation or divide current segment into 2 and
            % continue (Divide and Conquer)
            if max_val < line_acc_thres
                % Current line interpolation is accurate enough
                intvs = intv;
            else
                intv1 = obj.segmentation([intv(1) intv(1)+max_idx-1],depth+1);
                intv2 = obj.segmentation([intv(1)+max_idx-1 intv(2)],depth+1);
                intvs = [intv1 intv2];
            end
        
        end
        
        %% Merge line interpolated segments for valid circular approximation
        function obj = mergeSeg(obj,intvL)
            % Merge line intervals for valid circular approximation and
            % save data to "obj.init_segments"
            lb = 1; n = numel(intvL); cnt = 1;
            while true
                
                [~,err] = CircleFit(obj.LP(1,intvL(lb):intvL(n))',...
                                    obj.LP(2,intvL(lb):intvL(n))',...
                                    obj.w(intvL(lb):intvL(n))',...
                                    obj.cov(:,intvL(lb):intvL(n)),...
                                    0.99,false);
                if err.valid
                    disp(['Segment No. ',num2str(cnt),...
                          ' Idx ',num2str(intvL(lb)),' ~ ',num2str(intvL(n))])
                    res.bnds = [lb n];
                    obj.init_segments = [obj.init_segments res];
                    break;
                end

                max_fit_err = 0;
                % Error threshold value should be adequately large for safe
                % merging. If too large, there may not be appropriate 
                % intervals for good quality data, and if too small, short 
                % noisy intervals may be considered as a segment for this step.
                % Current sensor fusion model merge lanes almost perfectly,
                % so setting err_thres to about 1m will be fine.
                
                err_thres = 1; 
                while max_fit_err < err_thres
                    % Upper bound is found randomly between (lb+1,n) so that
                    % maximum circular fitting error is between values 3 and 5
                    ub = randint(lb+1,n);                    
                    [~,err] = CircleFit(obj.LP(1,intvL(lb):intvL(ub))',...
                                        obj.LP(2,intvL(lb):intvL(ub))',...
                                        obj.w(intvL(lb):intvL(ub))',...
                                        obj.cov(:,intvL(lb):intvL(ub)),...
                                        0.99,false);
                    max_fit_err = err.emax;
                end
    
                valid = err.valid;
    
                while ~valid
                    ub = ub - 1;
                    
                    if ub == lb
                        error('Need more line segmentation. Lower the line acc value')
                    end
                    [res,err] = CircleFit(obj.LP(1,intvL(lb):intvL(ub))',...
                                          obj.LP(2,intvL(lb):intvL(ub))',...
                                          obj.w(intvL(lb):intvL(ub))',...
                                          obj.cov(:,intvL(lb):intvL(ub)),...
                                          0.99,false);
                    valid = err.valid;
                    
                end
                disp(['Segment No. ',num2str(cnt),...
                      ' Idx ',num2str(intvL(lb)),' ~ ',num2str(intvL(ub))])
                res.bnds = [intvL(lb) intvL(ub)];
                obj.init_segments = [obj.init_segments res];
                cnt = cnt + 1;
                lb = ub;
                
                if lb == n
                    error('Need more line segmentation. Lower the line acc value')
                end
            end
            obj.num_seg = cnt;

        end

    end
end

%% ===========================Functions===========================

%% De-normalizing Constraints
function ER = InvMahalanobis(Xdiff, Cov)
    % Inverse Mahalanobis Distance for converting NLS problem to LS 
    n = size(Cov,1);
    SIG = eye(n)/chol(Cov);
    SIG = SIG';
    ER = SIG * Xdiff;
end

%% Get random integer with bnds
function val = randint(val1,val2)
    % find val from [val1, val2] randomly
    val = val1-1+randi(val2-val1+1);
end

%% Classify Lane Points based on angles
function assoc = classifyLP(LP,segments,init_flag)
    n = length(segments);
    assoc = {};
    % LP variable's 3rd row indicate absolute indices
    
    if init_flag
        for i=1:n
            bnds = segments{i}.bnds;
            res = struct();                        
            res.LP_idxs = bnds(1):bnds(2);
            res.LP_points = LP(:,bnds(1):bnds(2));      
            assoc = [assoc {res}];
        end    
    else
        LP_cpyd = LP;
        for i=1:n
            res = struct();
            if i~=n
                x = segments{i}.xc;
                y = segments{i}.yc;
                th_init = segments{i}.th_init;
                th_last = segments{i}.th_last;
                ths = atan2(LP_cpyd(2,:) - y, LP_cpyd(1,:) - x);
                ths_valididx = findInBetweenIdxs(ths,[th_init,th_last],i);
                
                LP_cpyd_valid = LP_cpyd(:,ths_valididx);
                LP_cpyd(:,ths_valididx) = [];
            else
                LP_cpyd_valid = LP_cpyd;
            end
                                
            res.LP_idxs = LP_cpyd_valid(3,:);
            res.LP_points = LP_cpyd_valid(1:2,:); 
            assoc = [assoc {res}];
        end
    end
end

%% Determine if angle is between boundary values
function idxs = findInBetweenIdxs(vals,bnds,seg_num)
    if bnds(1) > bnds(2)
        bnds = flip(bnds);
    end
    idxs = false(1,length(vals));
    
    for i=1:length(vals)
        if vals(i) >= bnds(1) && vals(i) <= bnds(2)
            idxs(i) = true;
        end
    end
    % Because of numerical error, first data point is excluded during the
    % data association process. Prevent first data point from be associated
    % to final segment(leftover data points)
    if seg_num == 1 
        idxs(1) = true;
    end
end
