function [res, err] = Circle2DFitV3(LP_l,LP_r,cov_l,cov_r,thres,intvs,plot_flag)
% Circle2DFitV3: Best fit Circle given data points and covariances
% Validity of circle fit is computed using the Mahalanobis distance
% Multiple Circular Regression given segmentation data indices
% Uses "Circle2DFitV2.m" for computing validity of each segments at initial
% stage. By comparing adjacent segment delL value, we can know whether
% curvature sign should be changing or not. This fact is reflected during
% the implementation of this script.
% To avoid having nonlinear equality constraints, different variables are
% used for other remaining segments excluding the first segment.
% 
% Implemented by Jin Hwan Jeon, 2022
% Email: jordan98@kaist.ac.kr
    LP_l = LP_l(1:2,:); LP_r = LP_r(1:2,:);
    % Creating Intervals and check validity for all segments separately
    n = length(intvs);
    if n <= 2
        error(['Wrong segmentation input, please enter at least 3 index for appropriate fitting...',...
               'For Segment-wise separate fitting, use Circle2DFitV2.m'])
    end
    m = n-1;

    seg_intvs = zeros(m,2);
    disp('<Step 1: Creating intervals and checking validity for all segments independently>')
    init_segments = {};
    for j=1:m
        disp('===============================')
        disp(['Checking validity of Segment ',num2str(j)])
        seg_intvs(j,1) = intvs(j);
        seg_intvs(j,2) = intvs(j+1);
        disp(['Data Index: ',num2str(seg_intvs(j,1)),'~',num2str(seg_intvs(j,2))])
        [res_,err_] = Circle2DFitV2(LP_l(:,seg_intvs(j,1):seg_intvs(j,2)),...
                                    LP_r(:,seg_intvs(j,1):seg_intvs(j,2)),...
                                    cov_l(:,seg_intvs(j,1):seg_intvs(j,2)),...
                                    cov_r(:,seg_intvs(j,1):seg_intvs(j,2)),...
                                    thres,false);
        res_.bnds = [seg_intvs(j,1), seg_intvs(j,2)];
        disp(['R: ',num2str(res_.R),' delL: ',num2str(res_.delL)])
        if j == 1
            init_delL = res_.delL;
            prev_delL = init_delL;
            res_.delL_status = 'init';
            res_.ang_status = 'init';
        else
            if res_.delL * init_delL < 0
                res_.delL_status = 'negative'; % delL sign should be different for this segment
            else
                res_.delL_status = 'positive'; % delL sign is kept equal for this segment
            end

            if res_.delL * prev_delL < 0
                res_.ang_status = 'negative'; % Angle of the first starting point should be pi - (previous segment last index angle)
            else
                res_.ang_status = 'positive'; % Angle of the first starting point should be (previous segment last index angle)
            end
            prev_delL = res_.delL;
        end
        init_segments = [init_segments {res_}];
        if err_.validity
            disp('Valid approximation')
        else
            disp('Invalid approximation')
        end
    end
    disp('===============================')
    disp('<Step 2: Multi-Arc Spline>')
    options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,...
                           'Display','iter-detailed','Algorithm','trust-region-reflective',...
                           'FunctionTolerance',1e-12,'MaxIterations',inf,'CheckGradients',false);

    % Compute Initial Values
    for j=1:m
        seg = init_segments{j};
        if j == 1
            X0_params = [seg.x; seg.y; seg.R; seg.delL];
            X0_th = seg.th';
            prev_R = seg.R;            
        else
            if strcmp(seg.ang_status,'positive')
                if seg.R > prev_R
                    X0_params = [X0_params; prev_R - seg.R]; 
                else
                    X0_params = [X0_params; seg.R - prev_R];
                end
                X0_th = [X0_th; seg.th(2:end)']; 
            else
                X0_params = [X0_params; prev_R + seg.R]; 
                X0_th = [X0_th; seg.th(2:end)']; 
            end
            
            % Exclude overlapping angle variables
            % Instead of radius, relative distance L is designed for
            % remaining segments to avoid nonlinear equality constraints.            
        end
    end
    X0 = [X0_params; X0_th];
    rel_base = seg_intvs(1,1); % Convert boundaries into relative indices
    % For example, intv = [400 500 600];
    % rel_idxs = [1 101 201];

    opt = lsqnonlin(@cost_func,X0,[],[],options);

    % Post Processing Optimization Results
    res = struct();
    res.segs = init_segments;

    opt_params = opt(1:4+m-1)';
    opt_th = opt(5+m-1:end)';
    
    res.opt_params = opt_params;
    res.opt_th = opt_th;

    for j=1:m
        lb = seg_intvs(j,1) - rel_base + 1; 
        ub = seg_intvs(j,2) - rel_base + 1;
        
        if j == 1           
            res.segs{j}.x = opt_params(1); res.segs{j}.y = opt_params(2); 
            res.segs{j}.R = opt_params(3); res.segs{j}.delL = opt_params(4);
            res.segs{j}.th = opt_th(lb:ub);
            prev_th = res.segs{j}.th(end);
        else
            % opt_x format
            % [L (th, excluded overlapping initial theta)]
            res.segs{j}.x = res.segs{j-1}.x + opt_params(4+j-1) * cos(opt_th(lb));
            res.segs{j}.y = res.segs{j-1}.y + opt_params(4+j-1) * sin(opt_th(lb));

            if strcmp(res.segs{j}.ang_status,'positive')                
                res.segs{j}.R = res.segs{j-1}.R - opt_params(4+j-1);
                res.segs{j}.th = opt_th(lb:ub); % overlapping initial angle is augmented
            else
                res.segs{j}.R = opt_params(4+j-1) - res.segs{j-1}.R;
                res.segs{j}.th = opt_th(lb:ub); % overlapping initial angle is augmented
                res.segs{j}.th(1) = opt_th(lb) + pi;
            end

            if strcmp(res.segs{j}.delL_status,'positive')
                res.segs{j}.delL = res.segs{1}.delL;
            else
                res.segs{j}.delL = -res.segs{1}.delL;
            end
        end

        [res.segs{j}.LP_l, res.segs{j}.LP_r] = get2DPoints(res.segs{j});
    end

    err = ComputeError();
    if plot_flag
        plotRes();
    end

    %% NLS Cost Function
    function [res, jac] = cost_func(x)       
        % Since sign of the curvature is not likely to change during
        % optimization, use init_segments to get coefficients used for 
        % converting constrained optimization problem to unconstrained
        % problem

        coeffs = zeros(m,2 + m-1); % Linear Coefficients for parameters
        
        
        % Pre-process to obtain coefficients for radius and beta(L_i's)
        for i=1:m 
            % coeffs format: [R_sign delL_sign (Li+1_signs)] 
            if i == 1                
                coeffs(i,1) = 1;
                coeffs(i,2) = -1;                
            else                
                seg_ = init_segments{i};
                
                if strcmp(seg_.delL_status,'positive')
                    coeffs(i,2) = -1;
                elseif strcmp(seg_.delL_status,'negative')
                    coeffs(i,2) = 1;
                end

                if strcmp(seg_.ang_status,'negative')
                    if i == 2
                        coeffs(i,i+1) = 1;
                        coeffs(i,1) = -1;
                    else
                        coeffs(i,[1,3:i]) = -coeffs(i-1,[1,3:i]);
                        coeffs(i,i+1) = 1;
                    end
                elseif strcmp(seg_.ang_status,'positive')
                    if i == 2
                        coeffs(i,i+1) = -1;
                        coeffs(i,1) = 1;
                    else
                        coeffs(i,[1,3:i]) = coeffs(i-1,[1,3:i]);
                        coeffs(i,i+1) = -1;
                    end
                end
            end
        end
        
        coeffs_l = coeffs;
        coeffs_r = coeffs;
        coeffs_r(:,2) = -coeffs(:,2);
        % Variable Clustering 
        x1 = x(1); y1 = x(2); R1 = x(3); delL = x(4);
        Ls = x(4+1:4+m-1)'; 
        th = x(5+m-1:end)'; % Total theta variables for data points, used seg_intvs to slice data
        
        R_params = [R1; delL; Ls'];
        res = []; jac = [];
        intvs_ = intvs(2:end-1);
        bnds = intvs_ - rel_base + 1;
        th_bnd = th(bnds); 
        
        % Computing Residual and Jacobian for each segment
        for i=1:m
            lb_ = seg_intvs(i,1) - rel_base + 1;
            ub_ = seg_intvs(i,2) - rel_base + 1;
            R_l = coeffs_l(i,:) * R_params;
            R_r = coeffs_r(i,:) * R_params;            

            if i == 1
                x_ = x1; y_ = y1;
                prev_x = x1;
                prev_y = y1;   
                th_ = th(lb_:ub_);
                
            else
                % th(lb_) angle variable for the previous segment
                x_ = prev_x + Ls(i-1) * cos(th(lb_));
                y_ = prev_y + Ls(i-1) * sin(th(lb_));
                prev_x = x_;
                prev_y = y_;
                                
                seg_ = init_segments{i};
                if strcmp(seg_.ang_status,'positive')
                    th_ = th(lb_:ub_);
                elseif strcmp(seg_.ang_status,'negative')
                    th_ = th(lb_:ub_);
                    th_(1) = pi + th_(1);
                end
            end
            
            % Creating Residual and blocks

            n_=  length(th_);
            
            blk_height = 2*2*n_;
            blk_width = length(x);
            resi = zeros(blk_height,1);
            
            I = []; J = []; V = [];
            
            for k=1:n_
                lp_pred_l = [x_ + R_l * cos(th_(k));
                             y_ + R_l * sin(th_(k))];
                lp_l = LP_l(:,seg_intvs(i,1) + k - 1);
                covl = reshape(cov_l(:,seg_intvs(i,1) + k - 1),2,2);
                resi(2*k-1:2*k) = InvMahalanobis(lp_pred_l - lp_l, covl);
                
                lp_pred_r = [x_ + R_r * cos(th_(k));
                             y_ + R_r * sin(th_(k))];
                lp_r = LP_r(:,seg_intvs(i,1) + k - 1);
                covr = reshape(cov_r(:,seg_intvs(i,1) + k - 1),2,2);
                resi(2*n_+2*k-1:2*n_+2*k) = InvMahalanobis(lp_pred_r - lp_r, covr);

                jacME_l = zeros(2,blk_width);
                jacME_r = zeros(2,blk_width);
                if i == 1
                    % x y elements
                    jacME_l(1,1) = 1; jacME_l(2,2) = 1;
                    % R delL elements
                    jacME_l(1,3:4) = coeffs_l(i,1:2) * cos(th_(k));
                    jacME_l(2,3:4) = coeffs_l(i,1:2) * sin(th_(k));
                    % th elements
                    jacME_l(1,4+m-1+lb_ + k - 1) = -R_l * sin(th_(k));
                    jacME_l(2,4+m-1+lb_ + k - 1) = R_l * cos(th_(k));

                    % x y elements
                    jacME_r(1,1) = 1; jacME_r(2,2) = 1;
                    % R delL elements
                    jacME_r(1,3:4) = coeffs_r(i,1:2) * cos(th_(k));
                    jacME_r(2,3:4) = coeffs_r(i,1:2) * sin(th_(k));
                    % th elements
                    jacME_r(1,4+m-1+lb_ + k - 1) = -R_r * sin(th_(k));
                    jacME_r(2,4+m-1+lb_ + k - 1) = R_r * cos(th_(k));
                else
                    % coeffs format: [R_sign delL_sign (Li+1_signs)]                                         
                    % x y elements
                    jacME_l(1,1) = 1; jacME_l(2,2) = 1;
                    % R delL elements
                    jacME_l(1,3:4) = coeffs_l(i,1:2) * cos(th_(k));
                    jacME_l(2,3:4) = coeffs_l(i,1:2) * sin(th_(k));
                    % Li's elements
                    jacME_l(1,4+1:4+i-1) = coeffs_l(i,2+1:2+i-1) * cos(th_(k)) + cos(th_bnd(1:i-1));
                    jacME_l(2,4+1:4+i-1) = coeffs_l(i,2+1:2+i-1) * sin(th_(k)) + sin(th_bnd(1:i-1));
                    % th elements
                    jacME_l(1,4+m-1+lb_ + k - 1) = -R_l * sin(th_(k));
                    jacME_l(2,4+m-1+lb_ + k - 1) = R_l * cos(th_(k));
                    % th elements -end         
                    
                    jacME_l(1,4+m-1+ bnds(1:i-1)) = jacME_l(1,4+m-1+ bnds(1:i-1)) - Ls(1:i-1).* sin(th_bnd(1:i-1));
                    jacME_l(2,4+m-1+ bnds(1:i-1)) = jacME_l(2,4+m-1+ bnds(1:i-1)) + Ls(1:i-1).* cos(th_bnd(1:i-1));

                    % x y elements
                    jacME_r(1,1) = 1; jacME_r(2,2) = 1;
                    % R delL elements
                    jacME_r(1,3:4) = coeffs_r(i,1:2) * cos(th_(k));
                    jacME_r(2,3:4) = coeffs_r(i,1:2) * sin(th_(k));
                    % Li's elements
                    jacME_r(1,4+1:4+i-1) = coeffs_r(i,2+1:2+i-1) * cos(th_(k)) + cos(th_bnd(1:i-1));
                    jacME_r(2,4+1:4+i-1) = coeffs_r(i,2+1:2+i-1) * sin(th_(k)) + sin(th_bnd(1:i-1));
                    % th elements
                    jacME_r(1,4+m-1+lb_ + k - 1) = -R_r * sin(th_(k));
                    jacME_r(2,4+m-1+lb_ + k - 1) = R_r * cos(th_(k));
                    % th elements -end                    
                    jacME_r(1,4+m-1+ bnds(1:i-1)) = jacME_r(1,4+m-1+ bnds(1:i-1)) - Ls(1:i-1).* sin(th_bnd(1:i-1));
                    jacME_r(2,4+m-1+ bnds(1:i-1)) = jacME_r(2,4+m-1+ bnds(1:i-1)) + Ls(1:i-1).* cos(th_bnd(1:i-1));
                end

                [I1,J1,V1] = sparseFormat(2*k-1:2*k,1:blk_width,InvMahalanobis(jacME_l,covl));
                [I2,J2,V2] = sparseFormat(2*n_+2*k-1:2*n_+2*k,1:blk_width,InvMahalanobis(jacME_r,covr));
                I = [I I1 I2]; J = [J J1 J2]; V = [V V1 V2];
                
            end
            MEsubBlock = sparse(I,J,V,blk_height,blk_width);
            res = [res; resi]; jac = [jac; MEsubBlock];            
        end
        
    end

    %% Compute Error
    function seg_errs = ComputeError()
        disp('Computing Error and Validty of current segmentation')
        seg_errs = {};
        
        for i=1:m
            disp(['[Segment ',num2str(i),' Analysis]'])
            seg_ = res.segs{i};
            n = length(seg_.th);

            err = struct();
            err.wfull_l = zeros(1,n); err.wfull_r = zeros(1,n);
            err.full_l = zeros(1,n); err.full_r = zeros(1,n);
            err.wthres = sqrt(chi2inv(thres,2));
            cnt_l = 0; cnt_r = 0;
            
            err.violated_l_idx = [];
            err.violated_r_idx = [];
            for k=1:n
                diff_l = seg_.LP_l(:,k) - LP_l(:,seg_intvs(i,1) + k - 1);
                diff_r = seg_.LP_r(:,k) - LP_r(:,seg_intvs(i,1) + k - 1);
                covl = reshape(cov_l(:,seg_intvs(i,1) + k - 1),2,2);
                covr = reshape(cov_r(:,seg_intvs(i,1) + k - 1),2,2);
                err.wfull_l(k) = sqrt(diff_l' / covl * diff_l);
                err.full_l(k) = sqrt(diff_l' * diff_l);
                err.wfull_r(k) = sqrt(diff_r' / covr * diff_r);
                err.full_r(k) = sqrt(diff_r' * diff_r);
                if err.wfull_l(k) > err.wthres
                    cnt_l = cnt_l + 1;
                    err.violated_l_idx = [err.violated_l_idx k];
                end
                if err.wfull_r(k) > err.wthres
                    cnt_r = cnt_r + 1;
                    err.violated_r_idx = [err.violated_r_idx k];
                end
            end

            if cnt_l > 2 || cnt_r > 2
                err.validity = false;
                disp('Current circular fitting is not valid, consider other index intervals for accurate fitting')
            else
                err.validity = true;
                disp("Current circular fitting is valid, extend index interval for more segment extension")
            end

            seg_errs = [seg_errs {err}];
        end
    end

    %% Plot Results
    function plotRes()
        figure(1); hold on; grid on; axis equal;
        p_data_l = plot(LP_l(1,intvs(1):intvs(end)),LP_l(2,intvs(1):intvs(end)),'r.');
        p_data_r = plot(LP_r(1,intvs(1):intvs(end)),LP_r(2,intvs(1):intvs(end)),'g.');
        for i=1:m
            seg_ = res.segs{i};
            seg_err = err{i};
            p_approx = plot(seg_.LP_l(1,:),seg_.LP_l(2,:),'k--');
            plot(seg_.LP_r(1,:),seg_.LP_r(2,:),'k--');
            
            % Control Points
            cts = [seg_.LP_l(:,1) seg_.LP_l(:,end) seg_.LP_r(:,1) seg_.LP_r(:,end)];
            p_ct = plot(cts(1,:),cts(2,:),'bp');
            
            % Threshold Violated Points
            vio = [seg_.LP_l(:,seg_err.violated_l_idx) seg_.LP_r(:,seg_err.violated_l_idx)];
            if ~isempty(vio)
                p_vio = plot(vio(1,:),vio(2,:),'cs');
                legend([p_data_l,p_data_r,p_approx,p_ct,p_vio],...
                'Left Lane Data','Right Lane Data', ...
                'Arc Spline','Control Points','Violated Data Points')
            else
                legend([p_data_l,p_data_r,p_approx,p_ct],...
                'Left Lane Data','Right Lane Data', ...
                'Arc Spline','Control Points')
            end
        end
        xlabel('Global X'); ylabel('Global Y'); title('Circular Fitting')
        
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
