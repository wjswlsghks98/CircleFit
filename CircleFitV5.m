classdef CircleFitV5 < handle
% [CircleFitV5: Nonlinear Least Squares based Multi-arc spline]
% Using curvature and arc_length as parameters
% To implement 1D weighting, only the lateral direction covariance is
% extracted for circular fitting
% Curve is continuously connected (2D and heading angle continuity is guaranteed)
% Implemented by Jin Hwan Jeon, 2022
% email: jordan98@kaist.ac.kr
    properties
        LP % Data points
        cov % Data point covariance
        w % Data point 1D weight, computed at initialization process
        thres % Chi-Square test accuracy threshold: 0.9, 0.95, 0.99
        intv % Data point interval of interest
        num_seg % Number of target arc segments
        init_segments = {} % Initial arc segments
        opt = struct() % Optimized Results
        segments = {} % Data points divided into segments        
    end

    methods
        %% Constructor
        function obj = CircleFitV5(LP,cov,thres,intv)
            obj.LP = LP(1:2,:);
            obj.cov = cov;
            obj.thres = thres;
            obj.intv = intv;             
        end

        %% Initialize
        function obj = initialize(obj)
            % [Part 1: Extract lateral direction covariance value]
            % When fitting data into arc splines, the most important part
            % is lateral fitting error. Although directly extracting
            % lateral direction covariance cannot be considered as a
            % perfectly accurate weight, it will show adequate performance
            % for fitting data points of different "importance".
            % For the covariance data I use, lateral covariance value is
            % much smaller than the longitudinal one, so taking "min" will
            % easily lead to extraction of lateral variance information

            n = size(obj.cov,2);
            obj.w = zeros(1,n);
            for i=1:n
                cov_ = reshape(obj.cov(:,i),2,2);
                obj.w(i) = min(eig(cov_));
            end

            % [Part 2: Create initial value segments]
            % Incrementally perform circular fitting to find the minimum
            % number of segments and their initial values
            % Initial number of segments is determined by this process
            
            obj.findInitialSegments();                     

            % [Part 3: Create initial value for optimization]
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
        end
        
        %% Find Initial Segments
        function obj = findInitialSegments(obj)
            % Starting from the index 1, incrementally find maximum
            % possible length arc spline approximation
            disp('=========[Initial Segmentization]=========');
            lb = 1; 
            % upper bound for search index will be found randomly first and
            % will be decreased. This is to increase the searching
            % efficiency of the maximum length accurate arc spline
            % approximation. Upper bound index is set so that the maximum
            % fitting error is between 5 ~ 10 m
            n = size(obj.LP,2);
            
            cnt = 1;
            obj.opt.init_intvs = lb;
            base_delta = 1000;
            while true
                % Determine if current segment of interest is the last
                % segment
                
                [res_,err_] = CircleFit(obj.LP(1,lb:n)',obj.LP(2,lb:n)',obj.w(lb:n)',obj.cov(:,lb:n),obj.thres,false);
                if err_.valid
                    res_.bnds = [lb n];
                    obj.init_segments = [obj.init_segments {res_}];
                    disp(['Segment ',num2str(cnt),' Index ',num2str(lb),'~',num2str(n)])
                    obj.opt.init_intvs = [obj.opt.init_intvs n];
                    break;
                end
                
                max_err = 0;
                % If not last segment, find upper bound
                ub = lb;
                while max_err < 3 
                    
                    ub = ub + base_delta; 
                    
                    if ub > n
                        ub = n;
                        break;
                    end
                    [~,err_] = CircleFit(obj.LP(1,lb:ub)',obj.LP(2,lb:ub)',obj.w(lb:ub)',obj.cov(:,lb:ub),obj.thres,false);
                    max_err = err_.emax;

                end
                

                validity = false;
                while ~validity
                    ub = ub - 1;
                    [res_,err_] = CircleFit(obj.LP(1,lb:ub)',obj.LP(2,lb:ub)',obj.w(lb:ub)',obj.cov(:,lb:ub),obj.thres,false);
                    validity = err_.valid;                    
                end

                res_.bnds = [lb ub];
                obj.init_segments = [obj.init_segments {res_}];
                disp(['Segment ',num2str(cnt),' Index ',num2str(lb),'~',num2str(ub)])
                obj.opt.init_intvs = [obj.opt.init_intvs ub];
                lb = ub;
                cnt = cnt + 1;
            end
            disp('=========[Initial Segmentization Finished]=========');
            obj.num_seg = cnt;
        end

        %% Optimize
        function obj = optimize(obj)
            options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',false,...
                                   'Display','iter-detailed','Algorithm','trust-region-reflective',...
                                   'UseParallel',true,'MaxFunctionEvaluations',inf,...
                                   'StepTolerance',1e-8,'MaxIterations',inf,'CheckGradients',false);

            obj.opt.valid = false;
            fig_base = length(obj.opt.init_intvs)-1;            

            while ~obj.opt.valid
                % Need to set lower bounds for arc lengths
                lb1 = [-inf,-inf,-inf];
                lb2 = [-inf 20];
                lb = horzcat(lb1,repmat(lb2,1,obj.num_seg));
                
                disp(['[Batch Optimization Starts, number of segments: ',num2str(obj.num_seg),']'])
                [obj.opt.fullx,~,obj.opt.res,~,~,~,obj.opt.jac] = lsqnonlin(@obj.cost_func,obj.opt.X0,lb,[],options);
                
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
                obj.computeError();    
                obj.plotRes(floor(obj.num_seg - fig_base + 1));                
                
                if ~obj.opt.valid
                    obj.createNewX0();                    
                end
            end
        end

        %% Cost Function
        function [res,jac] = cost_func(obj,x0)
            obj.segments = {};
            % [Part 1: Precompute parameter jacobians and calculate circle centers]
            obj.precompJac(x0);            
            % [Part 2: Divide data points using initial segment data]
            obj.classifyLP(x0);                  
            % [Part 3: Compute full jacobian and residual]
            % Part 3-1: Normal Measurement Jacobian using chain rule
            obj.CreateMEBlock();
            % Part 3-2: Initial point and Last point Measurement Jacobian
            obj.CreateMEBlock2(x0);

            obj.Merge();
            res = obj.opt.res;
            jac = obj.opt.jac;                                 
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

        %% Classify Lane Points based on angles
        function obj = classifyLP(obj,x0)
            n = length(obj.segments);
            LP_cpyd = [obj.LP(:,obj.intv(1):obj.intv(end)); 
                       obj.intv(1):obj.intv(end)];
            obj.opt.x0 = x0;
            if obj.num_seg ~= obj.opt.init_num_seg
                for i=1:n
                    if i~=n
                        x = obj.segments{i}.xc;
                        y = obj.segments{i}.yc;
                        th_init = obj.segments{i}.th_init;
                        th_last = obj.segments{i}.th_last;
                        ths = atan2(LP_cpyd(2,:) - y, LP_cpyd(1,:) - x);
                        ths_valididx = findInBetweenIdxs(ths,[th_init,th_last]);
                        
                        LP_cpyd_valid = LP_cpyd(:,ths_valididx);
                        LP_cpyd(:,ths_valididx) = [];
                    else
                        LP_cpyd_valid = LP_cpyd;
                    end
                    
                    obj.segments{i}.num = i;
                    obj.segments{i}.LP_idxs = LP_cpyd_valid(3,:);
                    obj.segments{i}.LP_points = LP_cpyd_valid(1:2,:); 
                end
            elseif norm(obj.opt.initX0 - x0) < 1e-4 % Numerical Error
                for i=1:n
                    bnds = obj.init_segments{i}.bnds;
                    obj.segments{i}.num = i;
                    if i == 1
                        obj.segments{i}.LP_idxs = bnds(1):bnds(2);
                        obj.segments{i}.LP_points = obj.LP(:,bnds(1):bnds(2));
                    else
                        obj.segments{i}.LP_idxs = (bnds(1)+1):bnds(2);
                        obj.segments{i}.LP_points = obj.LP(:,(bnds(1)+1):bnds(2));
                    end
                end
            else
                for i=1:n
                    if i~=n
                        x = obj.segments{i}.xc;
                        y = obj.segments{i}.yc;
                        th_init = obj.segments{i}.th_init;
                        th_last = obj.segments{i}.th_last;
                        ths = atan2(LP_cpyd(2,:) - y, LP_cpyd(1,:) - x);
                        ths_valididx = findInBetweenIdxs(ths,[th_init,th_last]);
                        
                        LP_cpyd_valid = LP_cpyd(:,ths_valididx);
                        LP_cpyd(:,ths_valididx) = [];
                    else
                        LP_cpyd_valid = LP_cpyd;
                    end
                    
                    obj.segments{i}.num = i;
                    obj.segments{i}.LP_idxs = LP_cpyd_valid(3,:);
                    obj.segments{i}.LP_points = LP_cpyd_valid(1:2,:); 
                end
            end
        end
        
        %% Create Normal Measurement Jacobian
        function obj = CreateMEBlock(obj)
            m = length(obj.segments);
            
            obj.opt.ME_block = zeros(obj.intv(end)-obj.intv(1)+1,3+2*obj.num_seg);
            obj.opt.ME_vec = zeros(obj.intv(end)-obj.intv(1)+1,1);
            cnt = 0;
            for i=1:m
                n = length(obj.segments{i}.LP_points);

                xc = obj.segments{i}.xc;
                yc = obj.segments{i}.yc;
                kappa = obj.segments{i}.kappa;
%                 if n <= 1
%                     disp([i n])
%                 end
                for j = cnt+1:cnt+n
                    LP_idx = obj.segments{i}.LP_idxs(j-cnt);
                    cov_ = obj.w(LP_idx);
                    xk = obj.segments{i}.LP_points(1,j-cnt);
                    yk = obj.segments{i}.LP_points(2,j-cnt);
                    d = sqrt((xc - xk)^2 + (yc - yk)^2);
                    resi = d - abs(1/kappa);
                    obj.opt.ME_vec(j) = InvMahalanobis(resi,cov_);
                    
                    if kappa >= 0
                        row_jac = [(xc - xk)/d, (yc - yk)/d, 1/kappa^2];
                    else
                        row_jac = [(xc - xk)/d, (yc - yk)/d, -1/kappa^2];
                    end

                    obj.opt.ME_block(j,:) = InvMahalanobis(row_jac * obj.segments{i}.Cjac,cov_);
                end

                cnt = cnt + n;
            end            
            obj.opt.ME_block = sparse(obj.opt.ME_block);            
        end

        %% Create Intial, Last Point Measurement Jacobian
        function obj = CreateMEBlock2(obj,x0)
            obj.opt.ME_block2 = zeros(4,3+2*obj.num_seg);
            obj.opt.ME_vec2 = zeros(4,1);
            % Initial Point
            cov1 = reshape(obj.cov(:,obj.intv(1)),2,2);

            obj.opt.ME_vec2(1:2) = InvMahalanobis(x0(1:2) - obj.LP(:,obj.intv(1)),cov1);
            obj.opt.ME_block2(1:2,1:2) = InvMahalanobis(eye(2),cov1);

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
            % compute jacobian
            base_jac = obj.segments{end}.Cjac(1:2,:);            
            
            add_jac = zeros(2,3+2*obj.num_seg);
            add_jac(:,3) = [1/kappas(end) * cos(alpha);
                            1/kappas(end) * sin(alpha)];
            for i=1:obj.num_seg
                if i~=obj.num_seg
                    add_jac(:,3+2*i-1) = [Ls(i)/kappas(end) * cos(alpha);
                                          Ls(i)/kappas(end) * sin(alpha)];
                    add_jac(:,3+2*i) = [kappas(i)/kappas(end) * cos(alpha);
                                        kappas(i)/kappas(end) * sin(alpha)];
                else
                    add_jac(:,3+2*i-1) = [Ls(i)/kappas(end) * cos(alpha) - 1/kappas(end)^2 * sin(alpha);
                                          Ls(i)/kappas(end) * sin(alpha) + 1/kappas(end)^2 * cos(alpha)];
                    add_jac(:,3+2*i) = [cos(alpha);
                                        sin(alpha)];
                end
            end
            
            obj.opt.ME_block2(3:4,:) = InvMahalanobis(base_jac+add_jac,cov2);            
        end

        %% Merge
        function obj = Merge(obj)
            obj.opt.res = vertcat(obj.opt.ME_vec,...
                                  obj.opt.ME_vec2);
            obj.opt.jac = vertcat(obj.opt.ME_block,...
                                  obj.opt.ME_block2);
        end
        
        %% Compute Error for every segment
        function obj = computeError(obj)
            m = length(obj.opt.segments);
            obj.opt.valid = true;
            obj.opt.err = {};

            disp('Error Analysis')
            for i=1:m                
                err_ = struct();
                n = length(obj.opt.segments{i}.LP_idxs);
                err_.full = zeros(1,n);
                
                err_.wfull = zeros(1,n);
                err_.wthres = sqrt(chi2inv(obj.thres,2));
                err_.invalid_idxs = [];
                
                xc = obj.opt.segments{i}.xc;
                yc = obj.opt.segments{i}.yc;
                R = abs(1/obj.opt.segments{i}.kappa);
                
                cnt = 0;
                for j=1:n
                    lp_idx = obj.opt.segments{i}.LP_idxs(j);
                    lp = obj.opt.segments{i}.LP_points(:,j);
                    th = atan2(lp(2)-yc,lp(1)-xc);
                    lp_pred = [xc + R * cos(th);
                               yc + R * sin(th)];
                    diff = lp_pred - lp;
                    cov_ = reshape(obj.cov(:,lp_idx),2,2);
                    err_.full(j) = sqrt(diff' * diff);
                    err_.wfull(j) = sqrt(diff' / cov_ * diff);
                    
                    if err_.wfull(j) > err_.wthres
                        cnt = cnt + 1;
                        err_.invalid_idxs = [err_.invalid_idxs j];
                    end
                end
                err_.emax = max(err_.full);
                
                if cnt < 3
                    disp(['Segment ',num2str(i),': Valid, Maximum Error: ',num2str(err_.emax)])
                    err_.valid = true;
                else
                    
                    disp(['Segment ',num2str(i),': Invalid, Maximum Error: ',num2str(err_.emax)])
                    err_.valid = false;
                    obj.opt.valid = false;
                end

                obj.opt.err = [obj.opt.err {err_}];
            end
        end

        %% Plot results
        function plotRes(obj,n)
            m = length(obj.segments);
            figure(n); hold on; grid on; axis equal;
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

                p_approx = plot(LPs(1,:),LPs(2,:),'k--');
                p_data = plot(obj.opt.segments{i}.LP_points(1,:),...
                              obj.opt.segments{i}.LP_points(2,:),'r.');
                p_cts = plot(point(1,:),point(2,:),'bp');

                if ~isempty(obj.opt.err{i}.invalid_idxs)
                    idxs = obj.opt.err{i}.invalid_idxs;
                    p_vio = plot(obj.opt.segments{i}.LP_points(1,idxs),...
                                 obj.opt.segments{i}.LP_points(2,idxs),'cs');
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
        end
        
        %% Create New X0 (initial value) for introducing new segment
        function obj = createNewX0(obj)
            n = length(obj.opt.segments);
            % Initial Point and initial heading angle is set equal
            
            x0 = [];
            invalid_cnt = 0;
            for i=1:n
                if ~obj.opt.err{i}.valid % Current Segment is not valid
                    x0 = [x0; obj.opt.segments{i}.kappa; 1/2 * obj.opt.segments{i}.L;
                              obj.opt.segments{i}.kappa; 1/2 * obj.opt.segments{i}.L];    
                    invalid_cnt = invalid_cnt + 1;
                else
                    x0 = [x0; obj.opt.segments{i}.kappa; obj.opt.segments{i}.L];
                end
            end
            
            obj.opt.X0 = [obj.opt.fullx(1:3); x0];  
            obj.num_seg = obj.num_seg + invalid_cnt;
        end

    end
end

%% ===========================Functions===========================

%% Determine if angle is between boundary values
function idxs = findInBetweenIdxs(vals,bnds)
    if bnds(1) > bnds(2)
        bnds = flip(bnds);
    end
    idxs = false(1,length(vals));
    
    for i=1:length(vals)
        if vals(i) >= bnds(1) && vals(i) <= bnds(2)
            idxs(i) = true;
        end
    end
end

%% De-normalizing Constraints
function ER = InvMahalanobis(Xdiff, Cov)
    % Inverse Mahalanobis Distance for converting NLS problem to LS 
    n = size(Cov,1);
    SIG = eye(n)/chol(Cov);
    SIG = SIG';
    ER = SIG * Xdiff;
end
