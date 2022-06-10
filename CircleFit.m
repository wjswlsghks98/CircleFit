function [x, y, R, error] = CircleFit(X, Y, w, plot_flag)
%% Circular Fitting
% From "Best Fit Circles made Easy" by Han de Bruijin
% URL: http://www.alternatievewiskunde.nl/jaar2006/kromming.pdf
% Robert Israel's method implemented by MATLAB 
% Implemented by JinHwan Jeon, 2022

% [x, y, R, error] = CircleFit(X, Y, w)
% 
% Returns best fit circle given 2D data points and corresponding weights

%% Input (should be in column vector form)
% X: x coordinates of data points 
% Y: y coordinates of data points
% w: weighting factor for particular point
% plot_flag: boolean for plotting fit results

%% Output
% x: x coordinate of optimized circle's center point
% y: y coordinate of optimized circle's center point
% R: radius of optimized circle
% error: fitting error information

% For understanding the implementation below, please read the paper
% provided in the URL above
    mu_x = w' * X / sum(w);
    mu_y = w' * Y / sum(w);
    
    X_rel = X - mu_x;
    Y_rel = Y - mu_y;
    sig_xx = sum(w.* X_rel.^2)/ sum(w);
    sig_yy = sum(w.* Y_rel.^2)/ sum(w);
    sig_xy = sum(w.* X_rel.* Y_rel)/ sum(w);
    
    if sig_xx * sig_yy > sig_xy^2
        sig_xxx = sum(w.* X_rel.^3)/ sum(w);
        sig_xxy = sum(w.* X_rel.^2.* Y_rel)/ sum(w);
        sig_xyy = sum(w.* X_rel.* Y_rel.^2)/ sum(w);
        sig_yyy = sum(w.* Y_rel.^3)/ sum(w);
        A = [sig_xx sig_xy;
             sig_xy sig_yy];
        b = [sig_xxx + sig_xyy;
             sig_xxy + sig_yyy];
        
        opt = A \ b;
        
        x_rel = 1/2 * opt(1);
        y_rel = 1/2 * opt(2);
        C = sig_xx + sig_yy;

        R = sqrt(C + x_rel^2 + y_rel^2);
        x = mu_x + x_rel;
        y = mu_y + y_rel;

    else
        mu_x = w' * X / sum(w);
        mu_y = w' * Y / sum(w);
        sig_xx = sum(w.* X.^2)/ sum(w);
        sig_yy = sum(w.* Y.^2)/ sum(w);
        sig_xy = sum(w.* X.* Y)/ sum(w);
        sig_xxx = sum(w.* X.^3)/ sum(w);
        sig_xxy = sum(w.* X.^2.* Y)/ sum(w);
        sig_xyy = sum(w.* X.* Y.^2)/ sum(w);
        sig_yyy = sum(w.* Y.^3)/ sum(w);

        A = [sig_xx sig_xy mu_x;
             sig_xy sig_yy mu_y;
             mu_x   mu_y   1];
        
        b = [sig_xxx + sig_xyy;
             sig_xxy + sig_yyy;
             sig_xx + sig_yy];
    
        opt = A \ b;
        x = 1/2 * opt(1);
        y = 1/2 * opt(2);
        R = sqrt(opt(3) + x^2 + y^2);
    end
    if ~isreal(R)
        warning('Complex radius detected, adjust dataset')
    end
    %% Plot Results
    if plot_flag
        figure(1);
        p_data = plot(X, Y, 'ro'); hold on; grid on; axis equal; 

        th = 0:pi/100:2*pi;
        xc = R * cos(th) + x;
        yc = R * sin(th) + y;

        p_fit = plot(xc,yc, 'k--');

        xlabel('X'); ylabel('Y'); title('Circular Fitting');
        legend([p_data, p_fit],'Original dataset', 'Optimized Circle')
    end

    %% Compute Fitting Error
    D = sqrt((X - x).^2 + (Y - y).^2);
    error = struct();
    error.full = D - R; % Naive Error

    % Non-Weighted Error Computation 
    error.se = sum((error.full).^2);
    error.mse = error.se / length(D);
    error.rmse = sqrt(error.mse);

    % Weighted Error Computation
    error.wse = sum((error.full).^2./w);
    error.wmse = error.wse / length(D);
    error.wrmse = sqrt(error.wmse);

end
