function [x, y, R, delL, theta, error] = CircleFitV2(X_l, X_r, Y_l, Y_r, w_l, w_r, plot_flag)
%% Circular Fitting for 2 Parallel Lanes
% Upgraded version from CircleFit for 2 lanes
% Since 2 lanes are parallel, they should share the same circle
% Implemented by JinHwan Jeon, 2022

% [x, y, R, delL, error] = CircleFit(X_l, X_r, Y_l, Y_r, w_l, w_r)
% 
% Returns best fit circle given 2D data points and corresponding weights

%% Input (should be in column vector form)
% X_l, X_r: x coordinates of data points for left and right lanes 
% Y_l, Y_r: y coordinates of data points for left and right lanes
% w_l, w_r: weighting factor for particular point (left and right lanes)
% plot_flag: boolean for plotting fit results

%% Output
% x: x coordinate of optimized circle's center point
% y: y coordinate of optimized circle's center point
% R: radius of optimized circle, strictly positive
% delL: lane spacing value for 2 lanes, may vary in sign depending on the
% location of center point of the circle
% theta: effective center circle angle
% error: fitting error information    
    
    % Data shifting to minimize scaling error 
    X_mean = mean([X_l; X_r]);
    Y_mean = mean([Y_l; Y_r]);
    X_l = X_l - X_mean; X_r = X_r - X_mean;
    Y_l = Y_l - Y_mean; Y_r = Y_r - Y_mean;

    A = zeros(4,4);
    b = zeros(4,1);

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
    
    theta = struct();
    theta.lb = min([th_list_l th_list_r]); theta.ub = max([th_list_l th_list_r]);
    theta.ang = theta.ub - theta.lb; 

    th_lbidx = find(th == theta.lb);
    th_ubidx = find(th == theta.ub);

    %% Plot Results
    if plot_flag
        figure(1);
        p_data_l = plot(X_l, Y_l, 'r.'); hold on; grid on; axis equal; 
        p_data_r = plot(X_r, Y_r, 'g.');

        p_fit = plot(xf1(th_lbidx:th_ubidx),yf1(th_lbidx:th_ubidx), 'k--');
        p_cont = plot(xf1(th_lbidx),yf1(th_lbidx), 'bp');
        plot(xf1(th_ubidx),yf1(th_ubidx), 'bp');

        plot(xf2(th_lbidx:th_ubidx),yf2(th_lbidx:th_ubidx),'k--');
        plot(xf2(th_lbidx),yf2(th_lbidx), 'bp');
        plot(xf2(th_ubidx),yf2(th_ubidx), 'bp');

        xlabel('X'); ylabel('Y'); title('Circular Fitting');
        legend([p_data_l, p_data_r, p_fit, p_cont],...
               'Original dataset(Left Lane)', 'Original dataset(Right Lane)', ...
               'Optimized Arcs', 'Arc Control Points')
    end

    %% Compute Fitting Error
    D_l = sqrt((X_l - x).^2 + (Y_l - y).^2);
    D_r = sqrt((X_r - x).^2 + (Y_r - y).^2);
    error = struct();
    error.full_l = D_l - (R - delL); % Naive Error
    error.full_r = D_r - (R + delL);
    error.max_l = max(error.full_l);
    error.max_r = max(error.full_r);

    % Non-Weighted Error Computation 
    error.se_l = sum((error.full_l).^2);
    error.mse_l = error.se_l / length(D_l);
    error.rmse_l = sqrt(error.mse_l);

    error.se_r = sum((error.full_r).^2);
    error.mse_r = error.se_r / length(D_r);
    error.rmse_r = sqrt(error.mse_r);

    % Weighted Error Computation
    error.wse_l = sum((error.full_l).^2./w_l);
    error.wmse_l = error.wse_l / length(D_l);
    error.wrmse_l = sqrt(error.wmse_l);

    error.wse_r = sum((error.full_r).^2./w_r);
    error.wmse_r = error.wse_r / length(D_r);
    error.wrmse_r = sqrt(error.wmse_r);

end
