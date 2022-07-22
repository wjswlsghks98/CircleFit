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
            obj.LP = LP;
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
            delta = floor((obj.intv(end)-obj.intv(1))/(obj.num_seg-1));
            init_intvs = obj.intv(1):delta:obj.intv(end);
            
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
            % [Part 1: Divide data points using initial segment data]
            obj.segments = ClassifyLP(obj.LP(:,obj.intvs(1):obj.intvs(end)),...
                                      obj.init_segments,obj.intvs);
            % [Part 2: Precompute parameter jacobians]
            % [Part 3: Compute full jacobian and residual]
        end
    end
end

%% Functions
function segments = ClassifyLP(LP,init_segments,intvs)
% Classify data points into several segments
    segments = {};
    n = length(init_segments);
    LP_cpyd = [LP;intvs(1):intvs(end)];
    for i=1:n
        x = init_segments{i}.x;
        y = init_segments{i}.y;
        th_init = init_segments{i}.th_init;
        th_last = init_segments{i}.th_last;
        ths = atan2(LP_cpyd(2,:) - y, LP_cpyd(1,:) - x);
        
    end
end
