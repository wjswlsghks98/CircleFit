classdef Circle2DFitV4 < handle
% Circle2DFitV4: Best fit Circle given data points and covariances
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
        opt = struct()
        auto
        init_segments = {}
    end

    methods
        %% Constructor
        function obj = Circle2DFitV4(LP_l,LP_r,cov_l,cov_r,thres,intvs,plot_flag,manual)
            obj.LP_l = LP_l(1:2,:);
            obj.LP_r = LP_r(1:2,:);
            obj.cov_l = cov_l;
            obj.cov_r = cov_r;
            obj.thres = thres;
            obj.intvs = intvs;
            if length(obj.intvs) < 3
                error('Enter more than 2 Segments. For segment-wise fitting, use Circle2DFitV2.m')
            end

            obj.plot_flag = plot_flag;
            obj.auto = manual; 
            % If true, compute jacobian inside the cost function
            % If false, compute jacobian by the solver (numerical differentiation)
        end
        
        %% Optimization
        function obj = optimize(obj)
            % Phase 1
            [obj,cont_flag] = obj.optimizePh1();
            
            % Phase 2
            if cont_flag
                obj.optimizePh2();
            else
                obj.opt.resnorm = 1e20; 
                % If optimization stopped, return 10^20 as the fitting cost
                % Setting this value too small will result in wrong optimal
                % segmentation intervals.
            end
        end
        
        %% Optimization Phase 1: Segment-wise analysis
        function [obj,flag] = optimizePh1(obj)
            flag = true; 
            % If any segment is detected to be infeasible at the initial 
            % stage, stop optimization because there is very low 
            % possibility that the segmentation becomes feasible after
            % multi-segment optimization.

            disp('Phase 1: Segment-wise approximation validity analysis')
            n = length(obj.intvs);
            m = n-1;
            % Compute Segment Intervals and 
            % Perform Optimization for each segment intervals

            obj.opt.seg_intvs = zeros(m,2);
            obj.opt.coeffs_l = zeros(m,2+m-1); % Linear Coefficients for parameters
            % Pre-process to obtain coefficients for radius and beta(L_i's)

            for i=1:m
                obj.opt.seg_intvs(i,1) = obj.intvs(i);
                obj.opt.seg_intvs(i,2) = obj.intvs(i+1);
                disp('===============================')
                disp(['Segment ',num2str(i),' Index: ',num2str(obj.intvs(i)),'~',num2str(obj.intvs(i+1))])
                % For optimization, overlapping indices are excluded!

                [res_,err_] = Circle2DFitV2(obj.LP_l,obj.LP_r,obj.cov_l,obj.cov_r,...
                                            obj.thres,obj.opt.seg_intvs(i,:),false); 

                res_.bnds = [obj.opt.seg_intvs(i,1), obj.opt.seg_intvs(i,2)];
                disp(['R: ',num2str(res_.R),' delL: ',num2str(res_.delL)])                

                if i == 1
                    init_delL = res_.delL;
                    prev_delL = init_delL;
                    res_.delL_status = 'init';
                    res_.ang_status = 'init';
                    
                    % Initial Value
                    X0_params = [res_.x; res_.y; res_.R; res_.delL];
                    X0_th = res_.th';
                    prev_R = res_.R;

                    % Coefficients
                    obj.opt.coeffs_l(i,1) = 1; obj.opt.coeffs_l(i,2) = -1;
                else
                    if res_.delL * init_delL < 0
                        res_.delL_status = 'negative'; % delL sign should be different for this segment
                        obj.opt.coeffs_l(i,2) = 1;
                    else
                        res_.delL_status = 'positive'; % delL sign is kept equal for this segment
                        obj.opt.coeffs_l(i,2) = -1;
                    end
        
                    if res_.delL * prev_delL < 0
                        res_.ang_status = 'negative'; % Angle of the first starting point should be pi - (previous segment last index angle)
                        X0_params = [X0_params; prev_R + res_.R];
                        
                        if i == 2
                            obj.opt.coeffs_l(i,i+1) = 1;
                            obj.opt.coeffs_l(i,1) = -1;
                        else
                            obj.opt.coeffs_l(i,[1,3:i]) = -obj.opt.coeffs_l(i-1,[1,3:i]);
                            obj.opt.coeffs_l(i,i+1) = 1;
                        end

                    else
                        res_.ang_status = 'positive'; % Angle of the first starting point should be (previous segment last index angle)
                        if res_.R > prev_R
                            X0_params = [X0_params; prev_R - res_.R];
                        else
                            X0_params = [X0_params; res_R - prev_R];
                        end

                        if i == 2
                            obj.opt.coeffs_l(i,i+1) = -1;
                            obj.opt.coeffs_l(i,1) = 1;
                        else
                            obj.opt.coeffs_l(i,[1,3:i]) = obj.opt.coeffs_l(i-1,[1,3:i]);
                            obj.opt.coeffs_l(i,i+1) = -1;
                        end
                    end
                    X0_th = [X0_th; res_.th(2:end)'];
                    prev_delL = res_.delL;
                end

                obj.init_segments = [obj.init_segments {res_}];

                if err_.validity
                    disp('Valid approximation')
                else
                    disp('Invalid approximation')
                    disp('Optimziation will not continue because infeasibility is detected while computing intial segment parameters')
                    flag = false;
                end
            end
            
            obj.opt.X0 = [X0_params; X0_th];
            obj.opt.coeffs_r = obj.opt.coeffs_l;
            obj.opt.coeffs_r(:,2) = -obj.opt.coeffs_l(:,2);
            disp('===============================')            
        end

        %% Optimization Phase 2: Constrained Optimization
        function obj = optimizePh2(obj)
            disp('Phase 2: Multi-segment arc spline optimization')
            options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,...
                                   'Display','iter-detailed','Algorithm','trust-region-reflective',...
                                   'MaxIterations',inf,'CheckGradients',false);
            
            % Solve Nonlinear Least Squares Problem
            [obj.opt.fullx,obj.opt.resnorm,~,~,~,~,obj.opt.jac] = lsqnonlin(@obj.cost_func,obj.opt.X0,[],[],options);
            
            % Post-Processing optimized data 
            obj.opt.segs = obj.init_segments;
            obj.opt.info = obj.opt.jac' * obj.opt.jac;
            obj.opt.info_density = nnz(obj.opt.info)/numel(obj.opt.info);

            m = length(obj.init_segments);
            obj.opt.params = obj.opt.fullx(1:4+m-1)';
            obj.opt.th = obj.opt.fullx(5+m-1:end)';
            
            rel_base = obj.opt.seg_intvs(1,1);
            
            % Saving optimized data into usable form (separate segments)
            for i=1:m
                lb = obj.opt.seg_intvs(i,1) - rel_base + 1;
                ub = obj.opt.seg_intvs(i,2) - rel_base + 1;
                    
                if i == 1           
                    obj.opt.segs{i}.x = obj.opt.params(1); 
                    obj.opt.segs{i}.y = obj.opt.params(2); 
                    obj.opt.segs{i}.R = obj.opt.params(3); 
                    obj.opt.segs{i}.delL = obj.opt.params(4);
                    obj.opt.segs{i}.th = obj.opt.th(lb:ub);
                else
                    % opt_x format
                    % [L (th, excluded overlapping initial theta)]
                    obj.opt.segs{i}.x = obj.opt.segs{i-1}.x + obj.opt.params(4+i-1) * cos(obj.opt.th(lb));
                    obj.opt.segs{i}.y = obj.opt.segs{i-1}.y + obj.opt.params(4+i-1) * sin(obj.opt.th(lb));
        
                    if strcmp(obj.opt.segs{i}.ang_status,'positive')                
                        obj.opt.segs{i}.R = obj.opt.segs{i-1}.R - obj.opt.params(4+i-1);
                        obj.opt.segs{i}.th = obj.opt.th(lb:ub); % overlapping initial angle is augmented
                    else
                        obj.opt.segs{i}.R = obj.opt.params(4+i-1) - obj.opt.segs{i-1}.R;
                        obj.opt.segs{i}.th = obj.opt.th(lb:ub); % overlapping initial angle is augmented
                        obj.opt.segs{i}.th(1) = th_add(obj.opt.th(lb),pi);
                    end
        
                    if strcmp(obj.opt.segs{i}.delL_status,'positive')
                        obj.opt.segs{i}.delL = obj.opt.segs{1}.delL;
                    else
                        obj.opt.segs{i}.delL = -obj.opt.segs{1}.delL;
                    end
                end
        
                [obj.opt.segs{i}.LP_l, obj.opt.segs{i}.LP_r] = get2DPoints(obj.opt.segs{i});

            end

            obj.computeError();

            if obj.plot_flag
                obj.plotRes();
            end
        end
        
        %% Optimization Cost Function
        function [res, jac] = cost_func(obj,x0)

            obj.CreateMMBlock(x0);
            obj.CreateMEBlock(x0);
            obj.Merge();
            
            res = obj.opt.res;
            jac = obj.opt.jac;
        end
        
        %% Create Motion Model Block
        function obj = CreateMMBlock(obj,x0)
            m = length(obj.init_segments); % Number of segments
            x1 = x0(1); y1 = x0(2); R1 = x0(3); delL = x0(4);
            Ls = x0(4+1:4+m-1)'; th = x0(5+m-1:end)';

            R_params = [R1; delL; Ls'];
            rel_base = obj.opt.seg_intvs(1,1);
            cov = 0.03^2;
            
            res = []; jac = [];
            for i=1:m    
                seg = obj.init_segments{i};
                R_l0 = seg.R - seg.delL;
                R_r0 = seg.R + seg.delL;
                th0 = seg.th;
                R_l = obj.opt.coeffs_l(i,:) * R_params;
                R_r = obj.opt.coeffs_r(i,:) * R_params;
                
                lb = obj.opt.seg_intvs(i,1) - rel_base + 1;
                ub = obj.opt.seg_intvs(i,2) - rel_base + 1;
                if i == 1                    
                    x = x1; y = y1; prev_x = x; prev_y = y;
                    th_ = th(lb:ub);
                else                    
                    x = prev_x + Ls(i-1) * cos(th(lb-1));
                    y = prev_y + Ls(i-1) * sin(th(lb-1));
                    prev_x = x; prev_y = y;
                    
                    if strcmp(seg.ang_status,'positive')
                        th_ = th(lb:ub);
                    else
                        th_ = th(lb:ub);
                        th_(1) = th_add(th_(1),pi);
                        if abs(th_(2) - th_(1)) > pi/4
                            error('Angle Approximation Error!')
                        end
                    end                   
                end
                
                if length(th_)~=length(th0)
                    error('Number of angle elements are not matched')
                end

                n = length(th_);
                blk_height = 2*(n-1);
                blk_width = length(x0);
                res_ = zeros(blk_height,1);
                I = []; J = []; V = [];
                for j=1:n-1
                    dth = th_(j+1) - th_(j);
                    % Left Lane
                    dLz_l = R_l0 * (th0(j+1) - th0(j));
                    dLzpred_l = R_l * (th_(j+1) - th_(j));
                    
                    res_(j) = InvMahalanobis(dLzpred_l - dLz_l,cov);
                    [I1,J1,V1] = sparseFormat(j,3:4+i-1,InvMahalanobis(dth * obj.opt.coeffs_l(i,1:2+i-1),cov));
                    [I2,J2,V2] = sparseFormat(j,4+m-1+lb+j-1:4+m-1+lb+j,InvMahalanobis([-R_l R_l],cov));

                    % Right Lane
                    dLz_r = R_r0 * (th0(j+1) - th0(j));
                    dLzpred_r = R_r * (th_(j+1) - th_(j));

                    res_(n-1+j) = InvMahalanobis(dLzpred_r - dLz_r,cov);
                    [I3,J3,V3] = sparseFormat(n-1+j,3:4+i-1,InvMahalanobis(dth * obj.opt.coeffs_r(i,1:2+i-1),cov));
                    [I4,J4,V4] = sparseFormat(n-1+j,4+m-1+lb+j-1:4+m-1+lb+j,InvMahalanobis([-R_r R_r],cov));
                    
                    I = [I I1 I2 I3 I4];
                    J = [J J1 J2 J3 J4];
                    V = [V V1 V2 V3 V4];
                end
                MMsubBlock = sparse(I,J,V,blk_height,blk_width);
                res = [res; res_];
                jac = [jac; MMsubBlock];
            end
            obj.opt.MM_res = res;
            obj.opt.MM_block = jac;
        end

        %% Create Measurement Block
        function obj = CreateMEBlock(obj,x0)
            m = length(obj.init_segments); % Number of segments
            x1 = x0(1); y1 = x0(2); R1 = x0(3); delL = x0(4);
            Ls = x0(4+1:4+m-1)'; th = x0(5+m-1:end)';

            R_params = [R1; delL; Ls'];
            rel_base = obj.opt.seg_intvs(1,1);
            bnds = obj.intvs(2:end-1) - rel_base + 1;
            th_bnd = th(bnds);
            
            res = []; jac = [];

            % Testing: exclude overlapping theta 
            for i=1:m
                R_l = obj.opt.coeffs_l(i,:) * R_params;
                R_r = obj.opt.coeffs_r(i,:) * R_params;

                if i == 1
                    lb = obj.opt.seg_intvs(i,1) - rel_base + 1;
                    ub = obj.opt.seg_intvs(i,2) - rel_base + 1;
                    x = x1; y = y1; prev_x = x; prev_y = y;
                    th_ = th(lb:ub);
                else
                    lb = obj.opt.seg_intvs(i,1) - rel_base + 1 + 1;
                    ub = obj.opt.seg_intvs(i,2) - rel_base + 1;
                    x = prev_x + Ls(i-1) * cos(th(lb-1));
                    y = prev_y + Ls(i-1) * sin(th(lb-1));
                    prev_x = x; prev_y = y;
                    
                    th_ = th(lb:ub);
%                     seg = obj.init_segments{i};
%                     if strcmp(seg.ang_status,'positive')
%                         th_ = th(lb:ub);
%                     else
%                         th_ = th(lb:ub);
%                         th_(1) = th_add(th_(1),pi);
%                     end
                end

                n = length(th_);
                
                blk_height = 2*2*n;
                blk_width = length(x0);
                res_ = zeros(blk_height,1);
                I = []; J = []; V = [];

                for j=1:n
                    lp_pred_l = [x + R_l * cos(th_(j));
                                 y + R_l * sin(th_(j))];
                    lp_l = obj.LP_l(:,lb + rel_base - 1 + j-1);
                    covl = reshape(obj.cov_l(:,lb + rel_base - 1 + j-1),2,2);

                    res_(2*j-1:2*j) = InvMahalanobis(lp_pred_l - lp_l, covl);
    
                    lp_pred_r = [x + R_r * cos(th_(j));
                                 y + R_r * sin(th_(j))];
                    lp_r = obj.LP_r(:,lb + rel_base - 1 + j-1);
                    covr = reshape(obj.cov_r(:,lb + rel_base - 1 + j-1),2,2);

                    res_(2*n+2*j-1:2*n+2*j) = InvMahalanobis(lp_pred_r - lp_r, covr);

                    jacME_l1 = zeros(2,4+i-1); jacME_l2 = zeros(2,1);
                    jacME_r1 = zeros(2,4+i-1); jacME_r2 = zeros(2,1);

                    if i == 1
                        % x y elements
                        jacME_l1(1,1) = 1; jacME_l1(2,2) = 1;
                        % R delL elements
                        jacME_l1(1,3:4) = obj.opt.coeffs_l(i,1:2) * cos(th_(j));
                        jacME_l1(2,3:4) = obj.opt.coeffs_l(i,1:2) * sin(th_(j));
                        
                        % th elements
                        jacME_l2(1) = -R_l * sin(th_(j));
                        jacME_l2(2) = R_l * cos(th_(j));
    
                        % x y elements
                        jacME_r1(1,1) = 1; jacME_r1(2,2) = 1;
                        % R delL elements
                        jacME_r1(1,3:4) = obj.opt.coeffs_r(i,1:2) * cos(th_(j));
                        jacME_r1(2,3:4) = obj.opt.coeffs_r(i,1:2) * sin(th_(j));
                        % th elements
                        jacME_r2(1) = -R_r * sin(th_(j));
                        jacME_r2(2) = R_r * cos(th_(j));
    
                    else
                        % coeffs format: [R_sign delL_sign (Li+1_signs)]                                         
                        % x y elements
                        jacME_l1(1,1) = 1; jacME_l1(2,2) = 1;
                        % R delL elements
                        jacME_l1(1,3:4) = obj.opt.coeffs_l(i,1:2) * cos(th_(j));
                        jacME_l1(2,3:4) = obj.opt.coeffs_l(i,1:2) * sin(th_(j));
                        % Li's elements
                        jacME_l1(1,4+1:4+i-1) = obj.opt.coeffs_l(i,2+1:2+i-1) * cos(th_(j)) + cos(th_bnd(1:i-1));
                        jacME_l1(2,4+1:4+i-1) = obj.opt.coeffs_l(i,2+1:2+i-1) * sin(th_(j)) + sin(th_bnd(1:i-1));                        
                        % th elements
                        jacME_l2(1) = -R_l * sin(th_(j));
                        jacME_l2(2) = R_l * cos(th_(j));
                        % th elements -end
                        jacME_l3 = zeros(2,i-1); 
                        jacME_l3(1,:) = - Ls(1:i-1).* sin(th_bnd(1:i-1));
                        jacME_l3(2,:) = Ls(1:i-1).* cos(th_bnd(1:i-1));

                        % x y elements
                        jacME_r1(1,1) = 1; jacME_r1(2,2) = 1;
                        % R delL elements
                        jacME_r1(1,3:4) = obj.opt.coeffs_r(i,1:2) * cos(th_(j));
                        jacME_r1(2,3:4) = obj.opt.coeffs_r(i,1:2) * sin(th_(j));
                        % Li's elements
                        jacME_r1(1,4+1:4+i-1) = obj.opt.coeffs_r(i,2+1:2+i-1) * cos(th_(j)) + cos(th_bnd(1:i-1));
                        jacME_r1(2,4+1:4+i-1) = obj.opt.coeffs_r(i,2+1:2+i-1) * sin(th_(j)) + sin(th_bnd(1:i-1));                        
                        % th elements
                        jacME_r2(1) = -R_l * sin(th_(j));
                        jacME_r2(2) = R_l * cos(th_(j));
                        % th elements -end
                        jacME_r3 = zeros(2,i-1); 
                        jacME_r3(1,:) = - Ls(1:i-1).* sin(th_bnd(1:i-1));
                        jacME_r3(2,:) = Ls(1:i-1).* cos(th_bnd(1:i-1));                              
                    end

                    [I1,J1,V1] = sparseFormat(2*j-1:2*j,1:4+i-1,InvMahalanobis(jacME_l1,covl));
                    [I2,J2,V2] = sparseFormat(2*j-1:2*j,4+m-1+lb+j-1,InvMahalanobis(jacME_l2,covl));
                    if i~=1
                        [I3,J3,V3] = sparseFormat(2*j-1:2*j,4+m-1+ bnds(1:i-1),InvMahalanobis(jacME_l3,covl));                            
                    else
                        I3 = []; J3 = []; V3 = [];
                    end
                    [I4,J4,V4] = sparseFormat(2*n+2*j-1:2*n+2*j,1:4+i-1,InvMahalanobis(jacME_r1,covr));
                    [I5,J5,V5] = sparseFormat(2*n+2*j-1:2*n+2*j,4+m-1+lb+j-1,InvMahalanobis(jacME_r2,covr));
                    if i~=1
                        [I6,J6,V6] = sparseFormat(2*n+2*j-1:2*n+2*j,4+m-1+ bnds(1:i-1),InvMahalanobis(jacME_r3,covl));                            
                    else
                        I6 = []; J6 = []; V6 = [];
                    end
    
                    I = [I I1 I2 I3 I4 I5 I6]; 
                    J = [J J1 J2 J3 J4 J5 J6]; 
                    V = [V V1 V2 V3 V4 V5 V6];
                end
                MEsubBlock = sparse(I,J,V,blk_height,blk_width);
                res = [res;res_];
                jac = [jac; MEsubBlock];
            end

            obj.opt.ME_res = res;
            obj.opt.ME_block = jac;
        end
        
        %% Merging Residual and Jacobian Block
        function obj = Merge(obj)
            obj.opt.res = vertcat(obj.opt.MM_res,...
                                  obj.opt.ME_res);
            obj.opt.jac = vertcat(obj.opt.MM_block,...
                                  obj.opt.ME_block);
        end

        %% Compute Optimization Error
        function obj = computeError(obj)
            m = length(obj.init_segments);
            obj.opt.err = {};
            for i=1:m
                disp(['[Segment ',num2str(i),' Analysis]'])
                seg_ = obj.opt.segs{i};
                n = length(seg_.th);
    
                err_ = struct();
                err_.wfull_l = zeros(1,n); err_.wfull_r = zeros(1,n);
                err_.full_l = zeros(1,n); err_.full_r = zeros(1,n);
                err_.wthres = sqrt(chi2inv(obj.thres,2));
                cnt_l = 0; cnt_r = 0;
                
                err_.violated_l_idx = [];
                err_.violated_r_idx = [];
                for j=1:n
                    diff_l = seg_.LP_l(:,j) - obj.LP_l(:,obj.opt.seg_intvs(i,1) + j - 1);
                    diff_r = seg_.LP_r(:,j) - obj.LP_r(:,obj.opt.seg_intvs(i,1) + j - 1);
                    covl = reshape(obj.cov_l(:,obj.opt.seg_intvs(i,1) + j - 1),2,2);
                    covr = reshape(obj.cov_r(:,obj.opt.seg_intvs(i,1) + j - 1),2,2);
                    err_.wfull_l(j) = sqrt(diff_l' / covl * diff_l);
                    err_.full_l(j) = sqrt(diff_l' * diff_l);
                    err_.wfull_r(j) = sqrt(diff_r' / covr * diff_r);
                    err_.full_r(j) = sqrt(diff_r' * diff_r);
                    if err_.wfull_l(j) > err_.wthres
                        cnt_l = cnt_l + 1;
                        err_.violated_l_idx = [err_.violated_l_idx j];
                    end
                    if err_.wfull_r(j) > err_.wthres
                        cnt_r = cnt_r + 1;
                        err_.violated_r_idx = [err_.violated_r_idx j];
                    end
                end
    
                if cnt_l > 2 || cnt_r > 2
                    err_.validity = false;
                    disp('Current circular fitting is not valid, consider other index intervals for accurate fitting')
                else
                    err_.validity = true;
                    disp("Current circular fitting is valid, extend index interval for more segment extension")
                end
    
                obj.opt.err = [obj.opt.err {err_}];
            end
        end

        %% Plot Results
        function plotRes(obj)
            % Parametrization Visualization
            figure(1); hold on; grid on; axis equal;
            p_data_l = plot(obj.LP_l(1,obj.intvs(1):obj.intvs(end)),obj.LP_l(2,obj.intvs(1):obj.intvs(end)),'r.');
            plot(obj.LP_r(1,obj.intvs(1):obj.intvs(end)),obj.LP_r(2,obj.intvs(1):obj.intvs(end)),'r.');
            m = length(obj.init_segments);
            for i=1:m
                seg_ = obj.opt.segs{i};
                seg_err = obj.opt.err{i};
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
            n = length(obj.opt.err);
            figure(2); hold on; grid on;
            for i=1:n
                err_ = obj.opt.err{i};
                x = obj.opt.segs{i}.bnds;
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

%% Angle addition
function th = th_add(th1, th2)
    th = th1 + th2 - floor((th1 + th2)/(2*pi)) * 2 * pi;
end
