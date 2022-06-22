function [res, err] = LineFitV2(X_l, X_r, Y_l, Y_r, w_l, w_r, plot_flag)
%% Line Fitting for 2 Parallel Line Data given
% Implemented by JinHwan Jeon, 2022
    
    A = zeros(3,3);
    A(1,1) = sum(w_l.* X_l.^2) + sum(w_r.* X_r.^2);
    A(1,2) = sum(w_l.* X_l);
    A(1,3) = sum(w_r.* X_r);
    
    A(2,1) = sum(w_l.* X_l);
    A(2,2) = sum(w_l);
    
    A(3,1) = sum(w_r.* X_r);
    A(3,3) = sum(w_r);
    
    b = zeros(3,1);
    b(1) = sum(w_l.* X_l.* Y_l) + sum(w_r.* X_r.* Y_r);
    b(2) = sum(w_l.* Y_l);
    b(3) = sum(w_r.* Y_r);
    
    res = struct();
    res.opt = A \ b;
    
    m = res.opt(1); n_l = res.opt(2); n_r = res.opt(3);
    res.delL = 1/2 * abs(n_l - n_r)/sqrt(m^2+1);
    
    [~,min_idxL] = min(X_l);
    [~,max_idxL] = max(X_l);
    [~,min_idxR] = min(X_r);
    [~,max_idxR] = max(X_r);
    
    if X_l(1) < X_l(end)
        res.left_to_right = true;
        res.theta = atan2(m,1);
        
        approxx_init_l = (X_l(min_idxL) + m * Y_l(min_idxL) - m*n_l)/(m^2 + 1^2);
        approxx_init_r = (X_r(min_idxR) + m * Y_r(min_idxR) - m*n_r)/(m^2 + 1^2);

        approxx_last_l = (X_l(max_idxL) + m * Y_l(max_idxL) - m*n_l)/(m^2 + 1^2);
        approxx_last_r = (X_r(max_idxR) + m * Y_r(max_idxR) - m*n_r)/(m^2 + 1^2);
        
        x_lim_init_l = approxx_init_l + res.delL * sin(res.theta);
        x_lim_init_r = approxx_init_r - res.delL * sin(res.theta);
        
        x_lim_last_l = approxx_last_l + res.delL * sin(res.theta);
        x_lim_last_r = approxx_last_r - res.delL * sin(res.theta);
        
        x_lim_init = min(x_lim_init_l, x_lim_init_r);
        x_lim_last = max(x_lim_last_l, x_lim_last_r);

    else
        res.left_to_right = false;
        res.theta = atan2(-m,-1);

        approxx_init_l = (X_l(max_idxL) + m * Y_l(max_idxL) - m*n_l)/(m^2 + 1^2);
        approxx_init_r = (X_r(max_idxR) + m * Y_r(max_idxR) - m*n_r)/(m^2 + 1^2);

        approxx_last_l = (X_l(min_idxL) + m * Y_l(min_idxL) - m*n_l)/(m^2 + 1^2);
        approxx_last_r = (X_r(min_idxR) + m * Y_r(min_idxR) - m*n_r)/(m^2 + 1^2);
        
        x_lim_init_l = approxx_init_l + res.delL * sin(res.theta);
        x_lim_init_r = approxx_init_r - res.delL * sin(res.theta);
        
        x_lim_last_l = approxx_last_l + res.delL * sin(res.theta);
        x_lim_last_r = approxx_last_r - res.delL * sin(res.theta);
        
        x_lim_init = max(x_lim_init_l, x_lim_init_r);
        x_lim_last = min(x_lim_last_l, x_lim_last_r);
        
    end
    
    x_lim_init_l = x_lim_init - res.delL * sin(res.theta);
    x_lim_init_r = x_lim_init + res.delL * sin(res.theta);
    
    x_lim_last_l = x_lim_last - res.delL * sin(res.theta);
    x_lim_last_r = x_lim_last + res.delL * sin(res.theta);
    
    res.init_pointL = [x_lim_init_l; m*x_lim_init_l + n_l];
    res.init_pointR = [x_lim_init_r; m*x_lim_init_r + n_r];
    
    res.last_pointL = [x_lim_last_l; m*x_lim_last_l + n_l];
    res.last_pointR = [x_lim_last_r; m*x_lim_last_r + n_r];
    
    diff = res.init_pointL - res.last_pointL;
    res.L = sqrt(diff' * diff);
    % Compute Error
    err = struct();
    err.tot = sum(w_l.* (m*X_l + n_l - Y_l).^2) + sum(w_r.* (m*X_r + n_r - Y_r).^2);
    err.full_l = abs(m * X_l - Y_l + n_l)/sqrt(m^2 + 1^2);
    err.full_r = abs(m * X_r - Y_r + n_r)/sqrt(m^2 + 1^2);
    err.max_l = max(err.full_l);
    err.max_r = max(err.full_r);
    
    err.se_l = sum((err.full_l).^2);
    err.mse_l = err.se_l / length(err.full_l);
    err.rmse_l = sqrt(err.mse_l);

    err.se_r = sum((err.full_r).^2);
    err.mse_r = err.se_r / length(err.full_r);
    err.rmse_r = sqrt(err.mse_r);
    
    % Plot Results
    if plot_flag
        left_lane = [res.init_pointL res.last_pointL];
        right_lane = [res.init_pointR res.last_pointR];
        
        figure(2);
        plot(X_l,Y_l,'r.'); hold on; grid on; axis equal;
        plot(X_r,Y_r,'g.');
        plot(left_lane(1,:),left_lane(2,:),'k--');
        plot(left_lane(1,:),left_lane(2,:),'bp');
        plot(right_lane(1,:),right_lane(2,:),'k--');
        plot(right_lane(1,:),right_lane(2,:),'bp');
        
    end
    
    
end
