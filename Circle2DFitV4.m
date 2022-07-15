classdef Circle2DFitV4 < handle
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

    properties
        LP_l
        LP_r
        cov_l
        cov_r
        thres
        intvs 
        plot_flag
        res
        err
    end

    methods
        %% Constructor
        function obj = Circle2DFitV4(LP_l,LP_r,cov_l,cov_r,thres,intvs,plot_flag)
            obj.LP_l = LP_l(1:2,:);
            obj.LP_r = LP_r(1:2,:);
            obj.cov_l = cov_l;
            obj.cov_r = cov_r;
            obj.thres = thres;
            obj.intvs = intvs;
            obj.plot_flag = plot_flag;
        end
        
        %% Optimize
        function obj = optimize(obj)
            n = length(obj.intvs);
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
                seg_intvs(j,1) = obj.intvs(j);
                seg_intvs(j,2) = obj.intvs(j+1);
                disp(['Data Index: ',num2str(seg_intvs(j,1)),'~',num2str(seg_intvs(j,2))])
                [res_,err_] = Circle2DFitV2(obj.LP_l(:,seg_intvs(j,1):seg_intvs(j,2)),...
                                            obj.LP_r(:,seg_intvs(j,1):seg_intvs(j,2)),...
                                            obj.cov_l(:,seg_intvs(j,1):seg_intvs(j,2)),...
                                            obj.cov_r(:,seg_intvs(j,1):seg_intvs(j,2)),...
                                            obj.thres,false);
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
                                   'FunctionTolerance',1e-12,'MaxIterations',inf,'CheckGradients',true);
        
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
            
            coeffs = zeros(m,2 + m-1); % Linear Coefficients for parameters
            % Pre-process to obtain coefficients for radius and beta(L_i's)
            for j=1:m 
                % coeffs format: [R_sign delL_sign (Li+1_signs)] 
                if j == 1                
                    coeffs(j,1) = 1;
                    coeffs(j,2) = -1;                
                else                
                    seg_ = init_segments{j};
                    
                    if strcmp(seg_.delL_status,'positive')
                        coeffs(j,2) = -1;
                    elseif strcmp(seg_.delL_status,'negative')
                        coeffs(j,2) = 1;
                    end
        
                    if strcmp(seg_.ang_status,'negative')
                        if j == 2
                            coeffs(j,j+1) = 1;
                            coeffs(j,1) = -1;
                        else
                            coeffs(j,[1,3:j]) = -coeffs(j-1,[1,3:j]);
                            coeffs(j,j+1) = 1;
                        end
                    elseif strcmp(seg_.ang_status,'positive')
                        if j == 2
                            coeffs(j,j+1) = -1;
                            coeffs(j,1) = 1;
                        else
                            coeffs(j,[1,3:j]) = coeffs(j-1,[1,3:j]);
                            coeffs(j,j+1) = -1;
                        end
                    end
                end
            end
            
            coeffs_l = coeffs;
            coeffs_r = coeffs;
            coeffs_r(:,2) = -coeffs(:,2);
            fin_jac = [];
            opt = lsqnonlin(@cost_func,X0,[],[],options);
            
            % Post Processing Optimization Results
            obj.res = struct();
            obj.res.segs = init_segments;
            obj.res.jac = fin_jac;
            obj.res.info = fin_jac' * fin_jac;
            obj.res.info_density = nnz(obj.res.info)/numel(obj.res.info);
        
            opt_params = opt(1:4+m-1)';
            opt_th = opt(5+m-1:end)';
            
            obj.res.opt_params = opt_params;
            obj.res.opt_th = opt_th;
        
            for j=1:m
                lb = seg_intvs(j,1) - rel_base + 1; 
                ub = seg_intvs(j,2) - rel_base + 1;
                
                if j == 1           
                    obj.res.segs{j}.x = opt_params(1); obj.res.segs{j}.y = opt_params(2); 
                    obj.res.segs{j}.R = opt_params(3); obj.res.segs{j}.delL = opt_params(4);
                    obj.res.segs{j}.th = opt_th(lb:ub);
                else
                    % opt_x format
                    % [L (th, excluded overlapping initial theta)]
                    obj.res.segs{j}.x = obj.res.segs{j-1}.x + opt_params(4+j-1) * cos(opt_th(lb));
                    obj.res.segs{j}.y = obj.res.segs{j-1}.y + opt_params(4+j-1) * sin(opt_th(lb));
        
                    if strcmp(obj.res.segs{j}.ang_status,'positive')                
                        obj.res.segs{j}.R = obj.res.segs{j-1}.R - opt_params(4+j-1);
                        obj.res.segs{j}.th = opt_th(lb:ub); % overlapping initial angle is augmented
                    else
                        obj.res.segs{j}.R = opt_params(4+j-1) - obj.res.segs{j-1}.R;
                        obj.res.segs{j}.th = opt_th(lb:ub); % overlapping initial angle is augmented
                        obj.res.segs{j}.th(1) = opt_th(lb) + pi;
                    end
        
                    if strcmp(obj.res.segs{j}.delL_status,'positive')
                        obj.res.segs{j}.delL = obj.res.segs{1}.delL;
                    else
                        obj.res.segs{j}.delL = -obj.res.segs{1}.delL;
                    end
                end
        
                [obj.res.segs{j}.LP_l, obj.res.segs{j}.LP_r] = get2DPoints(obj.res.segs{j});
            end
        
            obj.err = ComputeError();
            if obj.plot_flag
                obj.plotRes();
            end
            
            figure(3); spy(obj.res.info);
            disp('===============================')
            disp(['Density of Sparse Information Matrix: ',num2str(obj.res.info_density)])

        
            %% NLS Cost Function
            function [res, jac] = cost_func(x)       
                % Since sign of the curvature is not likely to change during
                % optimization, use init_segments to get coefficients used for 
                % converting constrained optimization problem to unconstrained
                % problem
                
                % Variable Clustering 
                x1 = x(1); y1 = x(2); R1 = x(3); delL = x(4);
                Ls = x(4+1:4+m-1)'; 
                th = x(5+m-1:end)'; % Total theta variables for data points, used seg_intvs to slice data
                
                R_params = [R1; delL; Ls'];
                res = []; jac = [];
                intvs_ = obj.intvs(2:end-1);
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
                        lp_l = obj.LP_l(:,seg_intvs(i,1) + k - 1);
                        covl = reshape(obj.cov_l(:,seg_intvs(i,1) + k - 1),2,2);
                        resi(2*k-1:2*k) = InvMahalanobis(lp_pred_l - lp_l, covl);
                        
                        lp_pred_r = [x_ + R_r * cos(th_(k));
                                     y_ + R_r * sin(th_(k))];
                        lp_r = obj.LP_r(:,seg_intvs(i,1) + k - 1);
                        covr = reshape(obj.cov_r(:,seg_intvs(i,1) + k - 1),2,2);
                        resi(2*n_+2*k-1:2*n_+2*k) = InvMahalanobis(lp_pred_r - lp_r, covr);
        
                        jacME_l1 = zeros(2,4+i-1); jacME_l2 = zeros(2,1);
                        jacME_r1 = zeros(2,4+i-1); jacME_r2 = zeros(2,1);
        
                        if i == 1
                            % x y elements
                            jacME_l1(1,1) = 1; jacME_l1(2,2) = 1;
                            % R delL elements
                            jacME_l1(1,3:4) = coeffs_l(i,1:2) * cos(th_(k));
                            jacME_l1(2,3:4) = coeffs_l(i,1:2) * sin(th_(k));
                            
                            % th elements
                            jacME_l2(1) = -R_l * sin(th_(k));
                            jacME_l2(2) = R_l * cos(th_(k));
        
                            % x y elements
                            jacME_r1(1,1) = 1; jacME_r1(2,2) = 1;
                            % R delL elements
                            jacME_r1(1,3:4) = coeffs_r(i,1:2) * cos(th_(k));
                            jacME_r1(2,3:4) = coeffs_r(i,1:2) * sin(th_(k));
                            % th elements
                            jacME_r2(1) = -R_r * sin(th_(k));
                            jacME_r2(2) = R_r * cos(th_(k));
        
        
                        else
                            
                            % coeffs format: [R_sign delL_sign (Li+1_signs)]                                         
                            % x y elements
                            jacME_l1(1,1) = 1; jacME_l1(2,2) = 1;
                            % R delL elements
                            jacME_l1(1,3:4) = coeffs_l(i,1:2) * cos(th_(k));
                            jacME_l1(2,3:4) = coeffs_l(i,1:2) * sin(th_(k));
                            % Li's elements
                            jacME_l1(1,4+1:4+i-1) = coeffs_l(i,2+1:2+i-1) * cos(th_(k)) + cos(th_bnd(1:i-1));
                            jacME_l1(2,4+1:4+i-1) = coeffs_l(i,2+1:2+i-1) * sin(th_(k)) + sin(th_bnd(1:i-1));
                            % th elements
                            if k == 1
                                jacME_l2(1) = -R_l * sin(th_(k)) - Ls(i-1) * sin(th_bnd(i-1));
                                jacME_l2(2) = R_l * cos(th_(k)) + Ls(i-1) * cos(th_bnd(i-1));
                            else
                                jacME_l2(1) = -R_l * sin(th_(k));
                                jacME_l2(2) = R_l * cos(th_(k));
                            end
                            % th elements -end         
                            if k == 1
                                jacME_l3 = zeros(2,i-2); 
                                jacME_l3(1,:) = - Ls(1:i-2).* sin(th_bnd(1:i-2));
                                jacME_l3(2,:) = Ls(1:i-2).* cos(th_bnd(1:i-2));
                            else
                                jacME_l3 = zeros(2,i-1); 
                                jacME_l3(1,:) = - Ls(1:i-1).* sin(th_bnd(1:i-1));
                                jacME_l3(2,:) = Ls(1:i-1).* cos(th_bnd(1:i-1));
                            end
                            
                            jacME_r1(1,1) = 1; jacME_r1(2,2) = 1;
                            % R delL elements
                            jacME_r1(1,3:4) = coeffs_r(i,1:2) * cos(th_(k));
                            jacME_r1(2,3:4) = coeffs_r(i,1:2) * sin(th_(k));
                            % Li's elements
                            jacME_r1(1,4+1:4+i-1) = coeffs_r(i,2+1:2+i-1) * cos(th_(k)) + cos(th_bnd(1:i-1));
                            jacME_r1(2,4+1:4+i-1) = coeffs_r(i,2+1:2+i-1) * sin(th_(k)) + sin(th_bnd(1:i-1));
                            % th elements
                            if k == 1
                                jacME_r2(1) = -R_r * sin(th_(k)) - Ls(i-1) * sin(th_bnd(i-1));
                                jacME_r2(2) = R_r * cos(th_(k)) + Ls(i-1) * cos(th_bnd(i-1));
                            else
                                jacME_r2(1) = -R_r * sin(th_(k));
                                jacME_r2(2) = R_r * cos(th_(k));
                            end
                            % th elements -end         
                            if k == 1
                                jacME_r3 = zeros(2,i-2);
                                jacME_r3(1,:) = - Ls(1:i-2).* sin(th_bnd(1:i-2));
                                jacME_r3(2,:) = Ls(1:i-2).* cos(th_bnd(1:i-2));
                            else
                                jacME_r3 = zeros(2,i-1);
                                jacME_r3(1,:) = - Ls(1:i-1).* sin(th_bnd(1:i-1));
                                jacME_r3(2,:) = Ls(1:i-1).* cos(th_bnd(1:i-1));
                            end
                        end
                        
                        [I1,J1,V1] = sparseFormat(2*k-1:2*k,1:4+i-1,InvMahalanobis(jacME_l1,covl));
                        [I2,J2,V2] = sparseFormat(2*k-1:2*k,4+m-1+lb_+k-1,InvMahalanobis(jacME_l2,covl));
                        if i~=1
                            if k == 1
                                [I3,J3,V3] = sparseFormat(2*k-1:2*k,4+m-1+ bnds(1:i-2),InvMahalanobis(jacME_l3,covl));
                            else
                                [I3,J3,V3] = sparseFormat(2*k-1:2*k,4+m-1+ bnds(1:i-1),InvMahalanobis(jacME_l3,covl));
                            end
                        else
                            I3 = []; J3 = []; V3 = [];
                        end
                        [I4,J4,V4] = sparseFormat(2*n_+2*k-1:2*n_+2*k,1:4+i-1,InvMahalanobis(jacME_r1,covr));
                        [I5,J5,V5] = sparseFormat(2*n_+2*k-1:2*n_+2*k,4+m-1+lb_+k-1,InvMahalanobis(jacME_r2,covr));
                        if i~=1
                            if k == 1
                                [I6,J6,V6] = sparseFormat(2*n_+2*k-1:2*n_+2*k,4+m-1+ bnds(1:i-2),InvMahalanobis(jacME_r3,covl));
                            else
                                [I6,J6,V6] = sparseFormat(2*n_+2*k-1:2*n_+2*k,4+m-1+ bnds(1:i-1),InvMahalanobis(jacME_r3,covl));
                            end
                        else
                            I6 = []; J6 = []; V6 = [];
                        end
        
                        I = [I I1 I2 I3 I4 I5 I6]; 
                        J = [J J1 J2 J3 J4 J5 J6]; 
                        V = [V V1 V2 V3 V4 V5 V6];
                        
                    end
                    MEsubBlock = sparse(I,J,V,blk_height,blk_width);
                    res = [res; resi]; jac = [jac; MEsubBlock];        
                end
                fin_jac = jac;
            end

            %% Compute Error
            function seg_errs = ComputeError()
                disp('Computing Error and Validty of current segmentation')
                seg_errs = {};
                
                for i=1:m
                    disp(['[Segment ',num2str(i),' Analysis]'])
                    seg_ = obj.res.segs{i};
                    n = length(seg_.th);
        
                    err__ = struct();
                    err__.wfull_l = zeros(1,n); err__.wfull_r = zeros(1,n);
                    err__.full_l = zeros(1,n); err__.full_r = zeros(1,n);
                    err__.wthres = sqrt(chi2inv(obj.thres,2));
                    cnt_l = 0; cnt_r = 0;
                    
                    err__.violated_l_idx = [];
                    err__.violated_r_idx = [];
                    for k=1:n
                        diff_l = seg_.LP_l(:,k) - obj.LP_l(:,seg_intvs(i,1) + k - 1);
                        diff_r = seg_.LP_r(:,k) - obj.LP_r(:,seg_intvs(i,1) + k - 1);
                        covl = reshape(obj.cov_l(:,seg_intvs(i,1) + k - 1),2,2);
                        covr = reshape(obj.cov_r(:,seg_intvs(i,1) + k - 1),2,2);
                        err__.wfull_l(k) = sqrt(diff_l' / covl * diff_l);
                        err__.full_l(k) = sqrt(diff_l' * diff_l);
                        err__.wfull_r(k) = sqrt(diff_r' / covr * diff_r);
                        err__.full_r(k) = sqrt(diff_r' * diff_r);
                        if err__.wfull_l(k) > err__.wthres
                            cnt_l = cnt_l + 1;
                            err__.violated_l_idx = [err__.violated_l_idx k];
                        end
                        if err__.wfull_r(k) > err__.wthres
                            cnt_r = cnt_r + 1;
                            err__.violated_r_idx = [err__.violated_r_idx k];
                        end
                    end
        
                    if cnt_l > 2 || cnt_r > 2
                        err__.validity = false;
                        disp('Current circular fitting is not valid, consider other index intervals for accurate fitting')
                    else
                        err__.validity = true;
                        disp("Current circular fitting is valid, extend index interval for more segment extension")
                    end
        
                    seg_errs = [seg_errs {err__}];
                end
            end
        end
        
        %% Plot Results
        function plotRes(obj)
            % Parametrization Visualization
            figure(1); hold on; grid on; axis equal;
            p_data_l = plot(obj.LP_l(1,obj.intvs(1):obj.intvs(end)),obj.LP_l(2,obj.intvs(1):obj.intvs(end)),'r.');
            plot(obj.LP_r(1,obj.intvs(1):obj.intvs(end)),obj.LP_r(2,obj.intvs(1):obj.intvs(end)),'r.');
            m = length(obj.res.segs);
            for i=1:m
                seg_ = obj.res.segs{i};
                seg_err = obj.err{i};
                p_approx = plot(seg_.LP_l(1,:),seg_.LP_l(2,:),'k--');
                plot(seg_.LP_r(1,:),seg_.LP_r(2,:),'k--');
                
                % Control Points
                cts = [seg_.LP_l(:,1) seg_.LP_l(:,end) seg_.LP_r(:,1) seg_.LP_r(:,end)];
                p_ct = plot(cts(1,:),cts(2,:),'bp');
                
                % Threshold Violated Points
                vio = [seg_.LP_l(:,seg_err.violated_l_idx) seg_.LP_r(:,seg_err.violated_r_idx)];
                vio_gt = [obj.LP_l(:,obj.intvs(i)+seg_err.violated_l_idx-1) obj.LP_r(:,obj.intvs(i)+seg_err.violated_r_idx-1)];
                if ~isempty(vio)
                    p_vio = plot(vio(1,:),vio(2,:),'cs');
                    plot(vio_gt(1,:),vio_gt(2,:),'ms')
                    legend([p_data_l,p_approx,p_ct,p_vio],...
                    'Lane Data',...
                    'Arc Spline','Control Points','Violated Data Points')
                else
                    legend([p_data_l,p_approx,p_ct],...
                    'Lane Data',...
                    'Arc Spline','Control Points')
                end
            end
            xlabel('Global X'); ylabel('Global Y'); title('Circular Fitting')
            
            % Parametrization Error Plot
            n = length(obj.err);
            figure(2); hold on; grid on;
            for i=1:n
                err_ = obj.err{i};
                x = obj.res.segs{i}.bnds;
                p_l = plot(x(1):x(2),err_.full_l,'r');
                p_r = plot(x(1):x(2),err_.full_r,'g');
            end
            xlabel('Index Number'); ylabel('Error(m)'); title('Fitting Error')
            legend([p_l,p_r],'Left Lane Error','Right Lane Error')
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
