function [res, err] = Circle2DFitV2(LP_l,LP_r,cov_l,cov_r,thres,plot_flag)
% Circle2DFit: Best fit Circle given data points and covariances
% Validity of circle fit is computed using the Mahalanobis distance (chi-square test)
% Angle variable reduced
% Implemented by Jin Hwan Jeon, 2022
% Email: jordan98@kaist.ac.kr

    options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,...
                           'Display','iter-detailed','Algorithm','trust-region-reflective',...
                           'FunctionTolerance',1e-12);
    LP_l = LP_l(1:2,:); LP_r = LP_r(1:2,:);
    X_l = LP_l(1,:); X_r = LP_r(1,:);
    Y_l = LP_l(2,:); Y_r = LP_r(2,:);
    n = length(X_l);
    % Compute Initial values
    [x0,y0,R0,~] = CircleFit(X_l',Y_l',ones(1,length(X_l))',false);
    X0 = [x0;y0;R0;0;zeros(1*length(X_l),1)];
    % a, b, R, delL, (th)
    opt = lsqnonlin(@cost_func,X0,[],[],options);
    res = struct();
    res.x = opt(1); res.y = opt(2); res.R = opt(3); res.delL = opt(4);
    res.th = opt(5:end)';
    [res.LP_l, res.LP_r] = get2DPoints(res);
    if res.th(1) > res.th(end)
        res.init_th = max(res.th);
        res.last_th = min(res.th);
    else
        res.init_th = min(res.th);
        res.last_th = max(res.th);
    end
    res.init_pointL = [res.x + (res.R - res.delL) * cos(res.init_th);
                       res.y + (res.R - res.delL) * sin(res.init_th)];
    res.last_pointL = [res.x + (res.R - res.delL) * cos(res.last_th);
                       res.y + (res.R - res.delL) * sin(res.last_th)];
    res.init_pointR = [res.x + (res.R + res.delL) * cos(res.init_th);
                       res.y + (res.R + res.delL) * sin(res.init_th)];
    res.last_pointR = [res.x + (res.R + res.delL) * cos(res.last_th);
                       res.y + (res.R + res.delL) * sin(res.last_th)];
    err = ComputeError();
    if plot_flag
        plotRes();
    end
    %% NLS Cost Function
    function [res, jac] = cost_func(x)
        a = x(1); b = x(2); R = x(3); delL = x(4);
        th = x(5:end);
        blk_height = 2*2*n;
        blk_width = 4 + n;
        res = zeros(blk_height,1);
        I = zeros(1,2*5*2*n); J = []; V = [];
        % Left Lane
        for i=1:n
            lp_pred = [a + (R - delL) * cos(th(i));
                       b + (R - delL) * sin(th(i))];
            lp = [X_l(i); Y_l(i)];
            cov = reshape(cov_l(:,i),2,2);
            [Jp, Jth] = JacME_l([a b R delL],th(i));
            res(2*i-1:2*i) = InvMahalanobis(lp_pred - lp,cov);
            [I1,J1,V1] = sparseFormat(2*i-1:2*i,1:4,InvMahalanobis(Jp,cov));
            [I2,J2,V2] = sparseFormat(2*i-1:2*i,4+i,InvMahalanobis(Jth,cov));
            I(10*i-9:10*i-2) = I1;J(10*i-9:10*i-2) = J1;V(10*i-9:10*i-2) = V1;
            I(10*i-1:10*i) = I2;J(10*i-1:10*i) = J2;V(10*i-1:10*i) = V2;
        end
        % Right Lane
        for i=1:n
            lp_pred = [a + (R + delL) * cos(th(i));
                       b + (R + delL) * sin(th(i))];
            lp = [X_r(i); Y_r(i)];
            cov = reshape(cov_r(:,i),2,2);
            [Jp, Jth] = JacME_r([a b R delL],th(i));
            res(2*n+2*i-1:2*n+2*i) = InvMahalanobis(lp_pred - lp,cov);
            [I1,J1,V1] = sparseFormat(2*n+2*i-1:2*n+2*i,1:4,InvMahalanobis(Jp,cov));
            [I2,J2,V2] = sparseFormat(2*n+2*i-1:2*n+2*i,4+i,InvMahalanobis(Jth,cov));
            I(10*n+10*i-9:10*n+10*i-2) = I1;J(10*n+10*i-9:10*n+10*i-2) = J1;V(10*n+10*i-9:10*n+10*i-2) = V1;
            I(10*n+10*i-1:10*n+10*i) = I2;J(10*n+10*i-1:10*n+10*i) = J2;V(10*n+10*i-1:10*n+10*i) = V2;
        end
        jac = sparse(I,J,V,blk_height,blk_width);
    end
    %% Compute Error
    function err = ComputeError()
        err = struct();
        err.wfull_l = zeros(1,n); err.wfull_r = zeros(1,n);
        err.full_l = zeros(1,n); err.full_r = zeros(1,n);
        err.wthres = sqrt(chi2inv(thres,2));
        cnt_l = 0; cnt_r = 0;
        for i=1:n
            diff_l = res.LP_l(:,i) - LP_l(:,i);
            diff_r = res.LP_r(:,i) - LP_r(:,i);
            covl = reshape(cov_l(:,i),2,2);
            covr = reshape(cov_r(:,i),2,2);
            err.wfull_l(i) = sqrt(diff_l' / covl * diff_l);
            err.full_l(i) = sqrt(diff_l' * diff_l);
            err.wfull_r(i) = sqrt(diff_r' / covr * diff_r);
            err.full_r(i) = sqrt(diff_r' * diff_r);
            if err.wfull_l(i) > err.wthres
                cnt_l = cnt_l + 1;
            end
            if err.wfull_r(i) > err.wthres
                cnt_r = cnt_r + 1;
            end
        end
        if cnt_l > 2 || cnt_r > 2
            err.validity = false;
            disp('Current circular fitting is not valid, consider other index intervals for accurate fitting')
        else
            err.validity = true;
            disp("Current circular fitting is valid, extend index interval for more segment extension")
        end
    end
    %% Plot Results
    function plotRes()
        figure(1); hold on; grid on; axis equal;
        p_data_l = plot(LP_l(1,:),LP_l(2,:),'r.');
        p_data_r = plot(LP_r(1,:),LP_r(2,:),'g.');
        p_approx = plot(res.LP_l(1,:),res.LP_l(2,:),'k--');
        plot(res.LP_r(1,:),res.LP_r(2,:),'k--');
        cts = [res.init_pointL res.last_pointL res.init_pointR res.last_pointR];
        p_ct = plot(cts(1,:),cts(2,:),'bp');
        xlabel('Global X'); ylabel('Global Y'); title('Circular Fitting')
        legend([p_data_l,p_data_r,p_approx,p_ct],...
                'Left Lane Data','Right Lane Data','Arc Spline','Control Points')
    end
end
%% Jacobians
function [Jp, Jth] = JacME_l(params,th)
    R = params(3); delL = params(4);
    Jp = [1 0 cos(th) -cos(th);
          0 1 sin(th) -sin(th)];
    Jth = (R-delL) * [-sin(th); cos(th)];
end

function [Jp, Jth] = JacME_r(params,th)
    R = params(3); delL = params(4);
    Jp = [1 0 cos(th) cos(th);
          0 1 sin(th) sin(th)];
    Jth = (R+delL) * [-sin(th); cos(th)];
end

%% De-normalizing Constraints
function ER = InvMahalanobis(Xdiff, Cov)
    % Inverse Mahalanobis Distance for converting NLS problem to LS
    n = size(Cov,1);
    SIG = eye(n)/chol(Cov);
    SIG = SIG';
    ER = SIG * Xdiff;
end

%% Return I, J, V for sparse matrix format
function [I, J, V] = sparseFormat(rows, cols, values)
    m = length(rows); n = length(cols);
    if ne(size(values), [m n])
        error('Value matrix format does not match row, column information')
    end
    I = zeros(1,m*n); J = zeros(1,m*n); V = zeros(1,m*n);
    for i=1:m
        for j=1:n
            I((i-1)*n+j) = rows(i);
            J((i-1)*n+j) = cols(j);
            V((i-1)*n+j) = values(i,j);
        end
    end
end

%% Get 2D Points
function [LP_l,LP_r] = get2DPoints(res)
    n = length(res.th);
    LP_l = zeros(2,n); LP_r = LP_l;
    for i=1:n
        LP_l(1,i) = res.x + (res.R - res.delL) * cos(res.th(i));
        LP_l(2,i) = res.y + (res.R - res.delL) * sin(res.th(i));
        LP_r(1,i) = res.x + (res.R + res.delL) * cos(res.th(i));
        LP_r(2,i) = res.y + (res.R + res.delL) * sin(res.th(i));
    end
end
