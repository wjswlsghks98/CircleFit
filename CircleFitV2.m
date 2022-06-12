function [x, y, R, delL, res, err] = CircleFitV2(varargin)
%% Circular Fitting for 2 Parallel Lanes
% Upgraded version from CircleFit for 2 lanes
% Since 2 lanes are parallel, they should share the same circle
% Implemented by JinHwan Jeon, 2022

% [x, y, R, delL, error] = CircleFit(X_l, X_r, Y_l, Y_r, w_l, w_r, plot_flag, init_flag, (lin_var: optional))
% 
% Returns best fit circle given 2D data points and corresponding weights

%% Input (should be in column vector form)
% X_l, X_r: x coordinates of data points for left and right lanes 
% Y_l, Y_r: y coordinates of data points for left and right lanes
% w_l, w_r: weighting factor for particular point (left and right lanes)
% plot_flag: boolean for plotting fit results
% init_flag: boolean for indicating current segment is intial or not
% lin_var: if not initial, perform constrained Least Squares using lin_var

%% Output
% x: x coordinate of optimized circle's center point
% y: y coordinate of optimized circle's center point
% R: radius of optimized circle, strictly positive
% delL: lane spacing value for 2 lanes, may vary in sign depending on the
% location of center point of the circle
% res: effective center circle angle + many more information
% error: fitting error information    
    
    %% Register Variables
    X_l = varargin{1};
    X_r = varargin{2};
    Y_l = varargin{3};
    Y_r = varargin{4};
    w_l = varargin{5};
    w_r = varargin{6};
    plot_flag = varargin{7};
    init_flag = varargin{8};

    if init_flag == true && nargin > 8
        error('To much input for initial arc segment fitting')
    elseif init_flag == false && nargin == 8
        error('Not enough input for arc segment fitting')
    end

    if ~init_flag
        init_points = varargin{9}; % Last Point information for previous segment
        prev_D = varargin{10}; % Previous Segment segment parameter D
        init_theta = varargin{11};
        % Used for determining current segment's lower/upper bound angle
        
        lp_l = init_points(:,1); lp_r = init_points(:,2);
    end

    %% Weighted Least Squares modeling for Circular Fitting
    % Data shifting to minimize scaling error 
    X_mean = mean([X_l; X_r]);
    Y_mean = mean([Y_l; Y_r]);
    X_l = X_l - X_mean; X_r = X_r - X_mean;
    Y_l = Y_l - Y_mean; Y_r = Y_r - Y_mean;
    
    if init_flag
        %% Initial Segment
        A = zeros(4,4); b = zeros(4,1);
    
        % Coefficient for "dQ/da = 0"
        A(1,1) = sum(w_l.* X_l.^2) + sum(w_r.* X_r.^2);
        A(1,2) = sum(w_l.* X_l.* Y_l) + sum(w_r.* X_r.* Y_r);
        A(1,3) = sum(w_l.* X_l) + sum(w_r.* X_r);
        A(1,4) = -sum(w_l.* X_l) + sum(w_r.* X_r);
    
        % Coefficient for "dQ/db = 0"
        A(2,1) = sum(w_l.* X_l.* Y_l) + sum(w_r.* X_r.* Y_r);
        A(2,2) = sum(w_l.* Y_l.^2) + sum(w_r.* Y_r.^2);
        A(2,3) = sum(w_l.* Y_l) + sum(w_r.* Y_r);
        A(2,4) = -sum(w_l.* Y_l) + sum(w_r.* Y_r);
    
        % Coefficient for "dQ/dC = 0"
        A(3,1) = sum(w_l.* X_l) + sum(w_r.* X_r);
        A(3,2) = sum(w_l.* Y_l) + sum(w_r.* Y_r);
        A(3,3) = sum(w_l) + sum(w_r);
        A(3,4) = -sum(w_l) + sum(w_r);
        
        % Coefficient for "dQ/dD = 0"
        A(4,1) = -sum(w_l.* X_l) + sum(w_r.* X_r);
        A(4,2) = -sum(w_l.* Y_l) + sum(w_r.* Y_r);
        A(4,3) = -sum(w_l) + sum(w_r);
        A(4,4) = sum(w_l) + sum(w_r);
        
        % b
        b(1) = sum(w_l.* X_l.^3) + sum(w_l.* X_l.* Y_l.^2) + ...
               sum(w_r.* X_r.^3) + sum(w_r.* X_r.* Y_r.^2);
    
        b(2) = sum(w_l.* X_l.^2.* Y_l) + sum(w_l.* Y_l.^3) + ...
               sum(w_r.* X_r.^2.* Y_r) + sum(w_r.* Y_r.^3);
    
        b(3) = sum(w_l.* X_l.^2) + sum(w_l.* Y_l.^2) + ...
               sum(w_r.* X_r.^2) + sum(w_r.* Y_r.^2);
    
        b(4) = -sum(w_l.* X_l.^2) - sum(w_l.* Y_l.^2) + ...
               sum(w_r.* X_r.^2) + sum(w_r.* Y_r.^2);
    
        opt = A \ b;
    
        x_rel = 1/2 * opt(1);
        y_rel = 1/2 * opt(2);
        C = opt(3); D = opt(4);
        
        % D > 0 ==> Circle is on the left lane region
        % D < 0 ==> Circle is on the right lane region
    
        RpdL = sqrt(x_rel^2 + y_rel^2 + C + D);
        
        x = x_rel + X_mean;
        y = y_rel + Y_mean;
    
        R = RpdL/2 + sqrt((RpdL/2)^2 - D/2);
        delL = RpdL/2 - sqrt((RpdL/2)^2 - D/2);
        
        if ~isreal(RpdL)
            warning('Incompatible results: adjust the sampling regions')
        end
    else
        %% Not Initial Segment
        x1 = lp_l(1) - X_mean; x2 = lp_r(1) - X_mean;
        y1 = lp_l(2) - Y_mean; y2 = lp_r(2) - Y_mean;
        m = (y2 - y1)/(x2 - x1);
        n = y1 - m * x1;

        A = zeros(3,3); b = zeros(3,1);

        % Coefficient for "dQ/da = 0"
        A(1,1) = sum(w_l.* (X_l + m*Y_l).^2) + sum(w_r.* (X_r + m*Y_r).^2);
        A(1,2) = sum(w_l.* (X_l + m*Y_l)) + sum(w_r.* (X_r + m*Y_r));
        A(1,3) = -sum(w_l.* (X_l + m*Y_l)) + sum(w_r.* (X_r + m*Y_r));

        % Coefficient for "dQ/dC = 0"
        A(2,1) = sum(w_l.* (X_l + m*Y_l)) + sum(w_r.* (X_r + m*Y_r));
        A(2,2) = sum(w_l) + sum(w_r);
        A(2,3) = -sum(w_l) + sum(w_r);

        % Coefficient for "dQ/dD = 0"
        A(3,1) = -sum(w_l.* (X_l + m*Y_l)) + sum(w_r.* (X_r + m*Y_r));
        A(3,2) = -sum(w_l) + sum(w_r);
        A(3,3) = sum(w_l) + sum(w_r);

        % b
        b(1) = sum(w_l.* (X_l + m*Y_l).* (X_l.^2 + Y_l.^2 - 2*n*Y_l)) + ...
               sum(w_r.* (X_r + m*Y_r).* (X_r.^2 + Y_r.^2 - 2*n*Y_r));
        
        b(2) = sum(w_l.* (X_l.^2 + Y_l.^2 - 2*n*Y_l)) + ...
               sum(w_r.* (X_r.^2 + Y_r.^2 - 2*n*Y_r));

        b(3) = -sum(w_l.* (X_l.^2 + Y_l.^2 - 2*n*Y_l)) + ...
               sum(w_r.* (X_r.^2 + Y_r.^2 - 2*n*Y_r));

        opt = A \ b;
        x_rel = 1/2 * opt(1);
        y_rel = m * x_rel + n;
        C = opt(2); D = opt(3);

        RpdL = sqrt(x_rel^2 + y_rel^2 + C + D);
        
        x = x_rel + X_mean;
        y = y_rel + Y_mean;
    
        R = RpdL/2 + sqrt((RpdL/2)^2 - D/2);
        delL = RpdL/2 - sqrt((RpdL/2)^2 - D/2);
        
        if ~isreal(RpdL)
            warning('Incompatible results: adjust the sampling regions')
        end
    end

    X_l = X_l + X_mean; X_r = X_r + X_mean;
    Y_l = Y_l + Y_mean; Y_r = Y_r + Y_mean;
    
    %% Compute Effective Angle 
    th = 0:pi/10000:2*pi;
    xf1 = (R - delL) * cos(th) + x;
    yf1 = (R - delL) * sin(th) + y;
    
    % Find approximate angle first and search around that angle
    d1 = (xf1 - X_l(1)).^2 + (yf1 - Y_l(1)).^2;
    d2 = (xf1 - X_l(end)).^2 + (yf1 - Y_l(end)).^2;

    [~,idx1] = min(d1);
    [~,idx2] = min(d2);
    mean_th = 1/2*(th(idx1) + th(idx2));
    delta = max(th(idx1),th(idx2)) - mean_th + 0.1;

    th = linspace(mean_th - delta, mean_th + delta, 1e5);
    
    xf1 = (R - delL) * cos(th) + x;
    yf1 = (R - delL) * sin(th) + y;

    xf2 = (R + delL) * cos(th) + x;
    yf2 = (R + delL) * sin(th) + y;

    th_list_l = zeros(1,length(X_l));
    th_list_r = th_list_l;

    for i=1:length(X_l)
        d_l = (xf1 - X_l(i)).^2 + (yf1 - Y_l(i)).^2;
        d_r = (xf2 - X_r(i)).^2 + (yf2 - Y_r(i)).^2;

        [~,idx_l] = min(d_l); [~,idx_r] = min(d_r);
        th_list_l(i) = th(idx_l); 
        th_list_r(i) = th(idx_r);
    end
    
    res = struct();
    res.D = D;
    % Approximated angle boundary calculation
    % If not initial segment, angle value for the intial data points will
    % be updated below
    res.th_lb = min([th_list_l th_list_r]); 
    res.th_ub = max([th_list_l th_list_r]);
    res.ang = res.th_ub - res.th_lb; 
    
    th_lbidx = find(th == res.th_lb);
    th_ubidx = find(th == res.th_ub);
    
    % Find the last point
    cnt_pointL_1 = [xf1(th_lbidx); yf1(th_lbidx)];
    cnt_pointL_2 = [xf1(th_ubidx); yf1(th_ubidx)];

    cnt_pointR_1 = [xf2(th_lbidx); yf2(th_lbidx)];
    cnt_pointR_2 = [xf2(th_ubidx); yf2(th_ubidx)];

    d1 = (X_l(end) - cnt_pointL_1(1))^2 + (Y_l(end) - cnt_pointL_1(2))^2;
    d2 = (X_l(end) - cnt_pointL_2(1))^2 + (Y_l(end) - cnt_pointL_2(2))^2;
    
    if d1 > d2
        res.last_pointL = cnt_pointL_2;
        res.last_pointR = cnt_pointR_2;
        
        res.next_theta = th(th_ubidx);
        
        % th(th_lbidx) is the initial point angle
        % res.th_lb should be changed if not initial segment
        
        if ~init_flag
            if D * prev_D > 0
                res.th_lb = init_theta;
            else
                res.th_lb = angleAdd(init_theta, pi);
            end
            res.ang = res.th_ub - res.th_lb;
        end
          
    else
        res.last_pointL = cnt_pointL_1;
        res.last_pointR = cnt_pointR_1;
        
        res.next_theta = th(th_lbidx);
        
        % th(th_ubidx) is the initial point angle
        % res.th_ub should be changed if not initial segment
        
        if ~init_flag
            if D * prev_D > 0
                res.th_ub = init_theta;
            else
                res.th_ub = angleAdd(init_theta, pi);
            end
            res.ang = res.th_ub - res.th_lb;
        end
    end
    
    new_th = linspace(res.th_lb, res.th_ub, 1e4);
    
    xf1 = (R - delL) * cos(new_th) + x;
    yf1 = (R - delL) * sin(new_th) + y;

    xf2 = (R + delL) * cos(new_th) + x;
    yf2 = (R + delL) * sin(new_th) + y;
    
    %% Plot Results
    if plot_flag
        figure(1);
        p_data_l = plot(X_l, Y_l, 'r.'); hold on; grid on; axis equal; 
        p_data_r = plot(X_r, Y_r, 'g.');

        p_fit = plot(xf1,yf1, 'k--');
        p_cont = plot(xf1(1),yf1(1), 'bp');
        plot(xf1(end),yf1(end), 'bp');

        plot(xf2,yf2,'k--');
        plot(xf2(1),yf2(1), 'bp');
        plot(xf2(end),yf2(end), 'bp');

        xlabel('X'); ylabel('Y'); title('Circular Fitting');
        legend([p_data_l, p_data_r, p_fit, p_cont],...
               'Original dataset(Left Lane)', 'Original dataset(Right Lane)', ...
               'Optimized Arcs', 'Arc Control Points')
    end

    %% Compute Fitting Error
    D_l = sqrt((X_l - x).^2 + (Y_l - y).^2);
    D_r = sqrt((X_r - x).^2 + (Y_r - y).^2);
    err = struct();
    
    % Non-Weighted Error Computation 
    err.full_l = abs(D_l - (R - delL)); 
    err.full_r = abs(D_r - (R + delL));
    err.max_l = max(err.full_l);
    err.max_r = max(err.full_r);

    err.se_l = sum((err.full_l).^2);
    err.mse_l = err.se_l / length(D_l);
    err.rmse_l = sqrt(err.mse_l);

    err.se_r = sum((err.full_r).^2);
    err.mse_r = err.se_r / length(D_r);
    err.rmse_r = sqrt(err.mse_r);

    % Weighted Error Computation
    err.wfull_l = err.full_l .* w_l.^(0.5);
    err.wfull_r = err.full_r .* w_r.^(0.5);
    err.wmax_l = max(err.wfull_l);
    err.wmax_r = max(err.wfull_r);

    err.wse_l = sum((err.full_l).^2.*w_l);
    err.wmse_l = err.wse_l / length(D_l);
    err.wrmse_l = sqrt(err.wmse_l);

    err.wse_r = sum((err.full_r).^2.*w_r);
    err.wmse_r = err.wse_r / length(D_r);
    err.wrmse_r = sqrt(err.wmse_r);

end

function sum = angleAdd(th1, th2)
    % Angle addition 
    n = floor((th1+th2)/(2*pi));
    sum = th1 + th2 - 2 * n * pi;
end
