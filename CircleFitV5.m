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
        function obj = CircleFitV5(LP,cov,thres,intv,num_seg)
            obj.LP = LP(1:2,:);
            obj.cov = cov;
            obj.thres = thres;
            obj.intv = intv; 
            obj.num_seg = num_seg;
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
            % Using the interval and number of segments, divide data points
            % and create initial segment data using "CircleFit.m" --> need
            % to modify
            init_intvs = floor(linspace(obj.intv(1),obj.intv(end),obj.num_seg+1));
%             disp(init_intvs)
            for i=1:obj.num_seg
                % Modify CircleFit output format
                % Need to calculate approximated arc_length and curvature
                [res,~] = CircleFit(obj.LP(1,init_intvs(i):init_intvs(i+1))',...
                                    obj.LP(2,init_intvs(i):init_intvs(i+1))',...
                                    obj.w(init_intvs(i):init_intvs(i+1))',...
                                    false);
                res.bnds = [init_intvs(i) init_intvs(i+1)];
                obj.init_segments = [obj.init_segments {res}];
            end

            % [Part 3: Create initial value for optimization]
            obj.opt.X0 = zeros(3+2*obj.num_seg,1);
            obj.opt.X0(1:2) = obj.LP(:,init_intvs(1));
            obj.opt.X0(3) = atan2(obj.LP(2,init_intvs(1)+1) - obj.LP(2,init_intvs(1)),...
                                  obj.LP(1,init_intvs(1)+1) - obj.LP(1,init_intvs(1)));
            for i=1:obj.num_seg
                obj.opt.X0(3+2*i-1:3+2*i) = [obj.init_segments{i}.kappa;
                                             obj.init_segments{i}.L];
            end
        end

        %% Optimize
        function obj = optimize(obj)
            options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,...
                                   'Display','iter-detailed','Algorithm','trust-region-reflective');
            obj.opt.fullx = lsqnonlin(@obj.cost_func,obj.opt.X0,[],[],options);
        end

        %% Cost Function
        function [res,jac] = cost_func(obj,x0)
            % [Part 1: Precompute parameter jacobians and calculate circle centers]
            obj.precompJac(x0);
            error('1');
            % [Part 2: Divide data points using initial segment data]
            obj.segments = ClassifyLP(obj.LP(:,obj.intv(1):obj.intv(end)),...
                                      obj.init_segments,obj.intv);
            
            
            % [Part 3: Compute full jacobian and residual]
        end
        
        %% Precompute Jacobians
        function obj = precompJac(obj,x0)
            tau = x0(3); kappa = x0(4); L = x0(5);
            xc = x0(1) - 1/kappa * sin(tau); yc = x0(2) + 1/kappa * cos(tau);
            jac_ = zeros(3,3+2*obj.num_seg);
            jac_(:,1:4) = [1, 0, -1/kappa * cos(tau), 1/kappa^2 * sin(tau);
                           0, 1, -1/kappa * sin(tau), -1/kappa^2 * cos(tau);
                           0, 0,                   0,                     1];
            seg = struct();
            seg.xc = xc; seg.yc = yc;
            seg.Cjac = jac_; % Chain Rule Based Pre-computed Jacobian
        
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
                    delta_jac(1,3+2*j) = kappaj/kappa * cos(tau);
                    delta_jac(2,3+2*j) = kappaj/kappa * sin(tau);
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
                    delta_jac(1,3+2*j) = delta_jac(1,3+2*j) - kappaj/kappa * cos(tau);
                    delta_jac(2,3+2*j) = delta_jac(2,3+2*j) - kappaj/kappa * sin(tau);
                end

                jac_ = jac_ + delta_jac;

                L = x0(3+2*i); % Li --> Li+1

                seg = struct();
                seg.xc = xc; seg.yc = yc;
                seg.Cjac = jac_; % Chain Rule Based Pre-computed Jacobian
            
                obj.segments = [obj.segments {seg}];
            end


        end
    end
end

%% ===========================Functions===========================
%% Classify Lane Points based on angles
function segments = ClassifyLP(LP,init_segments,intvs)
% Classify data points into several segments
    segments = {};
    n = length(init_segments);
    LP_cpyd = [LP;intvs(1):intvs(end)];
    for i=1:n
        if i~=n
            x = init_segments{i}.x;
            y = init_segments{i}.y;
            th_init = init_segments{i}.th_init;
            th_last = init_segments{i}.th_last;
            ths = atan2(LP_cpyd(2,:) - y, LP_cpyd(1,:) - x);
    
            ths_valididx = findInBetweenIdxs(ths,[th_init,th_last]);
            LP_cpyd_valid = LP_cpyd(:,ths_valididx);
            LP_cpyd(:,ths_valididx) = [];
        end

        seg = struct();
        seg.num = i;
        seg.LP_idxs = LP_cpyd_valid(3,:);
        seg.LP_points = LP_cpyd_valid(1:2,:);
        
        segments = [segments {seg}];
    end
end

%% Determine if angle is between boundary values
function idxs = findInBetweenIdxs(vals,bnds)
    if bnds(1) > bnds(2)
        bnds = flip(bnds);
    end
    
    idxs = false(1,length(vals));
    
    for i=1:length(vals)
        if vals(i) >= bnds(1) && vals(i) < bnds(2)
            idxs(i) = true;
        end
    end

end
