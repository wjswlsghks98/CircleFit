function [res, err] = CircleFit(X, Y, w, plot_flag)
%% Circular Fitting
% From "Best Fit Circles made Easy" by Han de Bruijin
% URL: http://www.alternatievewiskunde.nl/jaar2006/kromming.pdf
% Robert Israel's method implemented by MATLAB
% Implemented by JinHwan Jeon, 2022
% [res, err] = CircleFit(X, Y, w)
%
% Returns best fit circle given 2D data points and corresponding weights
%% Input (should be in column vector form)
% X: x coordinates of data points
% Y: y coordinates of data points
% w: weighting factor for particular point
% plot_flag: boolean for plotting fit results
%% Output [res, err]
% [res(struct)]: fitting result information
% res.x: x coordinate of optimized circle's center point
% res.y: y coordinate of optimized circle's center point
% res.R: radius of optimized circle
%
% [err(struct)]: fitting error information
% For understanding the implementation below, please read the paper
% provided in the URL above
    mu_x = w' * X / sum(w);
    mu_y = w' * Y / sum(w);
    X_rel = X - mu_x;
    Y_rel = Y - mu_y;
    sig_xx = sum(w.* X_rel.^2)/ sum(w);
    sig_yy = sum(w.* Y_rel.^2)/ sum(w);
    sig_xy = sum(w.* X_rel.* Y_rel)/ sum(w);

    res = struct();
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
    res.R = sqrt(C + x_rel^2 + y_rel^2);
    res.x = mu_x + x_rel;
    res.y = mu_y + y_rel;
    res.th_init = atan2(Y(1) - res.y, X(1) - res.x);
    res.th_last = atan2(Y(end) - res.y, X(end) - res.x);
    res.L = res.R * abs(res.th_last - res.th_init); % arc length

    if ~isreal(res.R)
        warning('Complex radius detected, adjust dataset')
    end
    %% Curvature Analysis
    % [Step 1: Determining initial data point heading angle]
    head_cnds = [res.th_init - pi/2, res.th_init + pi/2];
    diff = [X(1)-X(end);Y(1)-Y(end)] + (res.L/10) * [cos(head_cnds); sin(head_cnds)];
    d1 = diff(:,1)' * diff(:,1); d2 = diff(:,2)' * diff(:,2);
    
    if d1 < d2
        head = head_cnds(1);
    else
        head = head_cnds(2);
    end

    % [Step 2: Find sign of curvature]
    R_cnds = [res.R, -res.R];
    diff = [X(1)-res.x;Y(1)-res.y] + R_cnds.* [-sin(head);cos(head)];
    d1 = diff(:,1)' * diff(:,1); d2 = diff(:,2)' * diff(:,2);

    if d1 < d2
        res.kappa = 1/res.R;
    else
        res.kappa = -1/res.R;
    end

    %% Error Analysis
    err = struct();
%     D = sqrt((X - x).^2 + (Y - y).^2);
%     
%     err.full = D - R; % Naive Error
%     % Non-Weighted Error Computation
%     err.se = sum((err.full).^2);
%     err.mse = err.se / length(D);
%     err.rmse = sqrt(err.mse);
%     % Weighted Error Computation
%     err.wse = sum((err.full).^2./w);
%     err.wmse = err.wse / length(D);
%     err.wrmse = sqrt(err.wmse);

    %% Plot Results
    if plot_flag
        % Full Circle
        figure(1);
        p_data = plot(X, Y, 'ro'); hold on; grid on; axis equal;
        th = 0:pi/100:2*pi;
        xc = res.R * cos(th) + res.x;
        yc = res.R * sin(th) + res.y;
        p_fit = plot(xc,yc, 'k--');
        xlabel('X'); ylabel('Y'); title('Circular Fitting');
        legend([p_data, p_fit],'Original dataset', 'Optimized Circle')

        % Focused Plot
        figure(2);
        p_data = plot(X, Y, 'ro'); hold on; grid on; axis equal;
        th = linspace(res.th_last,res.th_init);
        xc = R * cos(th) + x;
        yc = R * sin(th) + y;
        p_fit = plot(xc,yc, 'k--');
        xlabel('X'); ylabel('Y'); title('Circular Fitting');
        legend([p_data, p_fit],'Original dataset', 'Optimized Circle')
    end
    
end
