classdef Optimizer_LinV4 < handle
% Optimizer module for Batch Optimization
% Lane Point is 1D lateral offset w.r.t vehicle position
% Since data association is not done properly, follow the following steps
% for stable optimization
% [Optimization Step]
% 1. INS propagation for initial value creation
% 2. Data Association (Finding maximum k)
% 3. NLS Optimization (Using lsqnonlin.m), converge with fixed data
% association from step 2
% 4. Compare residual norm with previous NLS optimization result 
% (If initial, previous optimized residual norm is set to inf).
% If difference in residual norm is smaller than certain threshold,
% stop. Else, increase counter and go to step 2. Stop also if counter goes
% over maximum iteration limit.
%
% Implemented by JinHwan Jeon, 2022
% contact: jordan98@kaist.ac.kr
    properties
        lane 
        output
        prev_num % Preview distance index 1~10
        lpc % Lane prediction constant 2~10
        dstart % Dataset start index
        dend % Dataset end index
        x0 % Optimization solver intial value
        lr = 1.44; % Distance from CoG(RT location) to rear axle
        % ls = 0.8; % Distance from CoG to IMU, Camera (Sensor Lever Arm)
        ins = struct(); % INS Mechanization Saver
        gnss = struct(); % GNSS Measurement Saver
        opt = struct(); % Optimzation Results Saver
        mode % Sensor Fusion Mode: 'full','partial'   
        plot_flag % Flag for automatic result plot
        partial % INS + GNSS + INS Fusion Results
    end
    methods
        %% Constructor
        function obj = Optimizer_LinV4(varargin)
            
            obj.lane = varargin{1};
            obj.output = varargin{2};
            obj.prev_num = varargin{3};
            obj.dstart = varargin{4};
            obj.dend = varargin{5};
            
            obj.mode = varargin{6};
            obj.plot_flag = varargin{7};
            
            if strcmp(obj.mode,'full')
                obj.partial = varargin{8};
            end

            dstart = obj.dstart;
            dend = obj.dend;            
            lane = obj.lane;
            
            if strcmp(obj.mode,'partial')
                obj.ins.dstart = dstart;
                obj.ins.dend = dend;
                obj.ins.iter = dstart;
                obj.ins.state = [lane.posx(dstart); lane.posy(dstart); lane.posz(dstart); 
                                 lane.vel.x(dstart); lane.vel.y(dstart); lane.vel.z(dstart);
                                 lane.eul.x(dstart); lane.eul.y(dstart); lane.eul.z(dstart);
                                 0; 0; 0;
                                 0; 0; 0];            
    
                obj.ins.state_saver = zeros(size(obj.ins.state,1),dend-dstart+1);
                obj.ins.state_saver(:,1) = obj.ins.state;                       
            
                
                
            elseif strcmp(obj.mode,'full')
                % For fast convergence, "partial" result is used.
                % This way, the optimization starts with much better
                % quality of data association
                obj.ins = obj.partial.ins;
            end
            obj.opt.lb = []; obj.opt.ub = [];
            obj.opt.options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,'Display','iter-detailed','Algorithm','trust-region-reflective');
            obj.opt.states = [];
            obj.opt.lml = [];
            obj.opt.lmr = [];
            
        end

        %% Initialize
        function obj = initialize(obj)
            % INS Mechanization for computing x0
            obj.INS();
            % Pre-process GNSS 
            obj.GNSS();
        end
        
        %% INS Mechanization
        function obj = INS(obj)
            if strcmp(obj.mode,'partial')
                disp('-INS Mechanization-')
                % INS Mechanization Loop
                while obj.ins.iter < obj.ins.dend
                    obj.ins.u = [obj.lane.xs.a(:,obj.ins.iter); obj.lane.xs.w(:,obj.ins.iter)];
                    obj.ins.dt = obj.lane.t(obj.ins.iter+1) - obj.lane.t(obj.ins.iter);
                    
                    obj.ins.state = StatePred(obj.ins.state, obj.ins.u, obj.ins.dt);               
                    obj.ins.iter = obj.ins.iter + 1;
                    
                    cnt = obj.ins.iter - obj.dstart + 1;
                    obj.ins.state_saver(:,cnt) = obj.ins.state;                
                end
            end                            
    
            % Merging vehicle states and lane points to create intial guess
            if strcmp(obj.mode,'full')
                % merge INS Mechanization results to create x0 (intial guess)
%                 obj.x0 = vertcat(1, reshape(obj.ins.state_saver,[],1), ...
%                                  reshape(obj.output.cam.l_inter(obj.dstart:obj.dend,1:obj.prev_num)',[],1), ...
%                                  reshape(obj.output.cam.r_inter(obj.dstart:obj.dend,1:obj.prev_num)',[],1));
            
                % 1 is wheel speed sensor scale factor initial value

                % Merge "partial" optimization results
                obj.x0 = vertcat(obj.partial.opt.x, ...
                                 reshape(obj.output.cam.l_inter(obj.dstart:obj.dend,1:obj.prev_num)',[],1), ...
                                 reshape(obj.output.cam.r_inter(obj.dstart:obj.dend,1:obj.prev_num)',[],1));
            elseif strcmp(obj.mode,'partial')
                obj.x0 = vertcat(1, reshape(obj.ins.state_saver,[],1));            
            else
                error('Options other than full, partial is not supported')
            end            
        end

        %% Pre-process GNSS measurements
        function obj = GNSS(obj)
            disp('-GNSS Pre-processing-')
            obj.gnss.lb = find(obj.output.ubloxgps.valid_laneidxs == obj.dstart);
            obj.gnss.ub = find(obj.output.ubloxgps.valid_laneidxs == obj.dend);
            hAcc_valid = obj.output.ubloxgps.hAcc(obj.output.ubloxgps.validity);
            if length(hAcc_valid) ~= length(obj.output.ubloxgps.valid_laneidxs)
                error('1')
            end
            vAcc_valid = obj.output.ubloxgps.vAcc(obj.output.ubloxgps.validity);
            x_valid = obj.output.ubloxgps.x(obj.output.ubloxgps.validity);
            y_valid = obj.output.ubloxgps.y(obj.output.ubloxgps.validity);
            z_valid = obj.output.ubloxgps.alt(obj.output.ubloxgps.validity) - obj.output.ubloxgps.alt_offset;

            if isempty(obj.gnss.lb) || isempty(obj.gnss.ub)
                if isempty(obj.gnss.lb) && ~isempty(obj.gnss.ub)
                    error('No valid GNSS measurement at start point only, select other data index for start point')
                elseif isempty(obj.gnss.lb) && isempty(obj.gnss.ub)
                    error('No valid GNSS measurement at both start and end points, select other data index for start and end points')
                elseif ~isempty(obj.gnss.lb) && isempty(obj.gnss.ub)
                    error('No valid GNSS measurement at end point only, select other data index for end point')
                end                
            end
            
%             if obj.output.ubloxgps.filtered_idxs(obj.gnss.lb) ~= obj.output.ubloxgps.scfull.trimmed_idxs(1) || ...
%                obj.output.ubloxgps.filtered_idxs(obj.gnss.ub) ~= obj.output.ubloxgps.scfull.trimmed_idxs(end)
%                 error('Trimmed GNSS indices need to be updated, default is 17301~22000')
%                 % To run optimization for different intervals, need to
%                 % correct "ubloxgps_proc.m" in Functions folder and save
%                 % updated results to "output.mat"
%             end

%             obj.gnss.cov = obj.output.ubloxgps.filtered_cov(:,obj.gnss.lb:obj.gnss.ub);
%             
%             obj.gnss.idxs = obj.output.ubloxgps.filtered_idxs(obj.gnss.lb:obj.gnss.ub);
%             % obj.gnss.iter_idxs = obj.gnss.idxs(obj.gnss.idxs <= obj.dstart);
%             
%             obj.gnss.meas = [obj.output.ubloxgps.filtered_x(obj.gnss.lb:obj.gnss.ub); 
%                              obj.output.ubloxgps.filtered_y(obj.gnss.lb:obj.gnss.ub);
%                              obj.output.ubloxgps.filtered_z(obj.gnss.lb:obj.gnss.ub)];   
            
            
            obj.gnss.idxs = obj.output.ubloxgps.valid_laneidxs(obj.gnss.lb:obj.gnss.ub);
            obj.gnss.stds = [hAcc_valid(obj.gnss.lb:obj.gnss.ub);
                             hAcc_valid(obj.gnss.lb:obj.gnss.ub);
                             vAcc_valid(obj.gnss.lb:obj.gnss.ub)];
            obj.gnss.meas = [x_valid(obj.gnss.lb:obj.gnss.ub);
                             y_valid(obj.gnss.lb:obj.gnss.ub);
                             z_valid(obj.gnss.lb:obj.gnss.ub)];
            figure(14); hold on; grid on; axis equal;
            for i=1:length(obj.gnss.idxs)
                gt = [obj.lane.posx(obj.gnss.idxs(i));
                      obj.lane.posy(obj.gnss.idxs(i))];
                ublx = obj.gnss.meas(1:2,i);

                pts = [gt ublx];
                plot(pts(1,:),pts(2,:));
                plot(gt(1),gt(2),'rx',ublx(1),ublx(2),'bx');
                
            end
            error('1')
        end

        %% Optimization via Solving Non-linear Least Squares
        function obj = optimize(obj)
            disp('-Batch Optimization-')
            disp(['-Current Mode: ',obj.mode,'-'])

            obj.opt.x = lsqnonlin(@obj.cost_func, obj.x0, obj.opt.lb, obj.opt.ub, obj.opt.options);
            
            [m,n] = size(obj.ins.state_saver);

            obj.opt.wsf = obj.opt.x(1);
            obj.opt.states = reshape(obj.opt.x(2:m*n+1),m,n);

            if strcmp(obj.mode,'full') 
                obj.opt.latL = reshape(obj.opt.x(1+m*n+1:1+m*n+n*obj.prev_num),obj.prev_num,n)';
                obj.opt.latR = reshape(obj.opt.x(1+m*n+n*obj.prev_num+1:end),obj.prev_num,n)';
                % need to convert 1D lateral distance data to 2D Lane Point
                % data
                obj.opt.lml = zeros(2*n,obj.prev_num);
                obj.opt.lmr = obj.opt.lml;
                for i=1:n
                    Xi = obj.opt.states(:,i);
                    for j=1:obj.prev_num
                        obj.opt.lml(2*i-1:2*i,j) = getLP(Xi,obj.opt.latL(i,j),j);
                        obj.opt.lmr(2*i-1:2*i,j) = getLP(Xi,obj.opt.latR(i,j),j);
                    end
                end
%                 obj.reorderLP();
            end
            
            obj.computeError();

            obj.opt.info_mat = obj.opt.jac' * obj.opt.jac;
            
            % Covariance Extraction using SuiteSparse Toolbox
            % Need to include for citation in paper
            
%             disp('-Recovering Covariance using SuiteSparse Toolbox-')
%             [obj.opt.cov, ~] = sparseinv(obj.opt.info_mat); 

            if obj.plot_flag
                obj.plotRes();
            end

%             if strcmp(obj.mode,'full') || strcmp(obj.mode,'2-phase')
%                 % Lane Point reordering
%                 obj.reorderLP();
%         
%             end
        end
        
        %% Optimization Cost Function
        function [res, jac] = cost_func(obj, x0)
            
            wsf = x0(1);
            m = size(obj.ins.state,1);
            n = size(obj.ins.state_saver,2);

            states = reshape(x0(2:m*n+1),m,n);
            
            if strcmp(obj.mode,'full') 
                latL = reshape(x0(m*n+1+1:m*n+1+n*obj.prev_num),obj.prev_num,n)';
                latR = reshape(x0(m*n+1+n*obj.prev_num+1:end),obj.prev_num,n)';
                lml = zeros(2*n,obj.prev_num);
                lmr = lml;
                for i=1:n
                    Xi = states(:,i);
                    for j=1:obj.prev_num
                        lml(2*i-1:2*i,j) = getLP(Xi,latL(i,j),j);
                        lmr(2*i-1:2*i,j) = getLP(Xi,latR(i,j),j);
                    end
                end
                
                obj.PreProcessing(states, lml, lmr); % Pre-Processing for every function iteration
            else
                latL = []; latR = [];
            end

            obj.CreatePrBlock(wsf, states, latL, latR);
            obj.CreateGNSSBlock(states); 
            obj.CreateMMBlock(states); 
            obj.CreateWSBlock(states, wsf); 
            obj.CreateMEBlock(states, latL, latR);
            
            obj.Merge();

            jac = obj.opt.jac;
            res = obj.opt.res;
           

            % Plot Optimization results every step
            figure(1);
            x = states(1,:); y = states(2,:);
            gnss_x = obj.gnss.meas(1,:); gnss_y = obj.gnss.meas(2,:);
            p_gnss = plot(gnss_x, gnss_y,'gx'); hold on; axis equal; grid on;
            p_opt2 = plot(x, y,'b.'); 
            p_opt1 = plot(x(obj.gnss.idxs-obj.dstart+1), y(obj.gnss.idxs-obj.dstart+1),'r.');
            
            if strcmp(obj.mode,'full') 
                                
                for i=1:obj.prev_num
                    opt_lml = reshape(lml(:,i),2,[]);
                    opt_lmr = reshape(lmr(:,i),2,[]);                        
                    p_left = plot(opt_lml(1,:), opt_lml(2,:));
                    p_right = plot(opt_lmr(1,:), opt_lmr(2,:));
                end
            end
            p_gt = plot(obj.lane.posx(obj.dstart:obj.dend), obj.lane.posy(obj.dstart:obj.dend),'k--');
            xlabel('Longitude(Deg)'); ylabel('Latitude(Deg)'); 
            
            
            if strcmp(obj.mode,'full')
                title('GNSS + INS + WSS + Lane Detection Sensor Fusion via MAP');
                legend([p_opt1 p_opt2 p_gt p_gnss p_left p_right],...
                       'Optimized Trajectory(w GNSS)','Optimized Trajectory(w/o GNSS)',...
                       'Ground Truth','GNSS Measurement(Ublox)',...
                       'Optimized Left Lane', 'Optimized Right Lane');
            
            elseif strcmp(obj.mode,'partial')
                title('GNSS + INS + WSS Sensor Fusion via MAP');
                legend([p_opt1, p_opt2, p_gt, p_gnss],...
                       'Optimized Trajectory(w GNSS)','Optimized Trajectory(w/o GNSS)',...
                       'Ground Truth','GNSS Measurement(Ublox)')
            end
            hold off;
        end
        
        %% Data Pre-Processing (Association)
        function obj = PreProcessing(obj, states, lml, lmr)
            disp('-Data Pre Processing-')
            
            [~, n] = size(states);
            
            obj.opt.assoc_idxsL = zeros(n,obj.prev_num);
            obj.opt.assoc_idxsR = zeros(n,obj.prev_num);
            for i=2:n
                
                for j=1:obj.prev_num
                    % Left Lane
                    k = 1;
                    
                    while true            
                        lp = lml(2*(i-k)-1:2*(i-k),j);                        
                        psi = states(9,i); % heading angle
                        delta = lp - states(1:2,i);                        
                        rot = [cos(psi) sin(psi)];
                        rel_x = rot * delta;                         
                        if rel_x < 0
                            k = k-1;
                            break;
                        else
                            k = k+1;
                            if i == k
                                k = k-1;
                                break;
                            end
                        end
                    end
                    obj.opt.assoc_idxsL(i,j) = k;
                    % Right Lane
                    k = 1;
                    
                    while true
                        lp = lmr(2*(i-k)-1:2*(i-k),j);                        
                        psi = states(9,i); % heading angle
                        delta = lp - states(1:2,i);                        
                        rot = [cos(psi) sin(psi)];
                        rel_x = rot * delta;
                        if rel_x < 0
                            k = k-1;
                            break;
                        else
                            k = k+1;
                            if i == k
                                k = k-1;
                                break;
                            end
                        end
                    end

                    obj.opt.assoc_idxsR(i,j) = k;
                end
            end
        end

        %% Prior Block and Residual
        function obj = CreatePrBlock(obj, wsf, states, latL, latR)
            % Initial Vehicle State
            prior_cov = diag([0.1^2 ... % WSS
                              0.01^2 0.01^2 0.01^2 ...
                              0.01^2 0.01^2 0.01^2 ...
                              0.01^2 0.01^2 0.01^2 ...
                              0.3^2 0.3^2 0.3^2 ...
                              0.1^2 0.1^2 0.1^2]);
    
            m = size(states,1);
            n = size(states,2);
            
            if strcmp(obj.mode,'full')
                blk_width = 1 + m*n + 2*n*obj.prev_num;
            elseif strcmp(obj.mode,'partial')
                blk_width = 1 + m*n;
            end
            
            % Vehicle State Prior
            Pr_vec_veh = zeros(m+1,1);
            Blk_s = sparse([],[],[],m+1,blk_width);
            Blk_s(:,1:m+1) = InvMahalanobis(eye(m+1),prior_cov);
            
            Pr_diff = vertcat(-1 + wsf, ...
                              -obj.ins.state_saver(:,1) + states(:,1));
            Pr_vec_veh(1:m+1) = InvMahalanobis(Pr_diff,prior_cov);
            
            if strcmp(obj.mode,'full') 
                % Lane Points
                Pr_vec_lane = zeros(2*n*obj.prev_num,1);
                I = zeros(1,2*n*obj.prev_num); J = I; V = I;
                

                for i=1:n
                    for j=1:obj.prev_num
                        cov_l = obj.output.cam.lstd_inter(obj.dstart+i-1,j)^2;
                        cov_r = obj.output.cam.rstd_inter(obj.dstart+i-1,j)^2;
                        
                        z_l = obj.output.cam.l_inter(obj.dstart+i-1,j);
                        z_r = obj.output.cam.r_inter(obj.dstart+i-1,j);
                        
                        zpred_l = latL(i,j);
                        zpred_r = latR(i,j);
                        Pr_vec_lane(obj.prev_num*(i-1)+j) = InvMahalanobis(zpred_l-z_l,cov_l);
                        Pr_vec_lane(obj.prev_num*(n+i-1)+j) = InvMahalanobis(zpred_r-z_r,cov_r);
                        I(obj.prev_num*(i-1)+j) = obj.prev_num*(i-1)+j;
                        I(obj.prev_num*(n+i-1)+j) = obj.prev_num*(n+i-1)+j;
                        J(obj.prev_num*(i-1)+j) = 1+m*n+obj.prev_num*(i-1)+j;
                        J(obj.prev_num*(n+i-1)+j) = 1+m*n+obj.prev_num*(n+i-1)+j;
                        V(obj.prev_num*(i-1)+j) = InvMahalanobis(1,cov_l);
                        V(obj.prev_num*(n+i-1)+j) = InvMahalanobis(1,cov_r);
                    end                                        
                end
                Blk_lp = sparse(I,J,V,n*2*obj.prev_num,blk_width);
                
                obj.opt.Pr_vec = vertcat(Pr_vec_veh, Pr_vec_lane);
                obj.opt.Pr_block = vertcat(Blk_s, Blk_lp);

            elseif strcmp(obj.mode,'partial')
                obj.opt.Pr_vec = Pr_vec_veh;
                obj.opt.Pr_block = Blk_s;
            end
        end
        
        %% GNSS Block and Residual
        function obj = CreateGNSSBlock(obj, states)
            Jx = horzcat(eye(3),zeros(3,12));
            [m,n] = size(states);
            mJ = size(Jx,1);

            if strcmp(obj.mode,'full')
                blk_width = m*n + 2*n*obj.prev_num;
            elseif strcmp(obj.mode,'partial')
                blk_width = m*n;
            end

            blk_height = size(obj.gnss.meas,2)*mJ;
            GNSS_vec = zeros(blk_height,1);
            
            tmp = m*3;
            d_gnss = size(obj.gnss.meas,2);
            
            I = zeros(1,tmp*d_gnss); J = I; V = I;

            for i=1:d_gnss
                state_idx = obj.gnss.idxs(i);
                rel_idx = state_idx - obj.dstart+1;
                
                gnss_cov = diag([obj.gnss.stds(1,i)^2, obj.gnss.stds(2,i)^2, obj.gnss.stds(3,i)^2]);
                             
                [I_s, J_s, V_s] = sparseFormat(mJ*i-mJ+1:mJ*i,m*(rel_idx-1)+1:m*rel_idx,InvMahalanobis(Jx,gnss_cov));
                I(tmp*i-(tmp-1):tmp*i) = I_s; J(tmp*i-(tmp-1):tmp*i) = J_s; V(tmp*i-(tmp-1):tmp*i) = V_s;
                
                % Residual
                Z_gnss = obj.gnss.meas(:,i);
                Z_pred = Jx * states(:,rel_idx);
                Z_diff = -Z_gnss + Z_pred;
                
                GNSS_vec(mJ*i-mJ+1:mJ*i) = InvMahalanobis(Z_diff, gnss_cov);
            end
            
            GNSS_block = sparse(I, J, V, blk_height, blk_width);
            dummy = sparse([],[],[], blk_height, 1);
            
            obj.opt.GNSS_vec = GNSS_vec;
            obj.opt.GNSS_block = horzcat(dummy, GNSS_block);

        end

        %% Motion Model Block and Residual
        function obj = CreateMMBlock(obj, states)
            [m,n] = size(states);
            
            MM_vec = zeros(m*(n-1),1);
            
            ab_std = 0.9e-5;
            wb_std = 2*1e-06;
            
            blk_height = m*(n-1);

            if strcmp(obj.mode,'full') 
                blk_width = m*n + 2*n*obj.prev_num;
            elseif strcmp(obj.mode,'partial')
                blk_width = m*n;
            end
            
            tmp = 2*m^2;
            I = zeros(1,tmp*(n-1)); J = I; V = I;
            for i=1:n-1
                
                u = [obj.lane.xs.a(:,obj.dstart+i-1); obj.lane.xs.w(:,obj.dstart+i-1)];
                dt = obj.lane.t(i+1) - obj.lane.t(i);

                ab_cov = 1/dt * ab_std^2;
                wb_cov = 1/dt * wb_std^2;
                
                MM_cov = diag([0.1^2 0.1^2 0.5^2 ...
                               0.05^2 0.05^2 0.1^2 ... 
                               0.005^2 0.005^2 0.005^2 ...
                               ab_cov ab_cov ab_cov ...
                               wb_cov wb_cov wb_cov]);

                X_curr = states(:,i); X_next = states(:,i+1);
                X_pred = StatePred(X_curr, u, dt);
                
                Jx = JacobianMM(X_curr, u, dt);
                   
                [I_s1, J_s1, V_s1] = sparseFormat(m*i-(m-1):m*i,m*i-(m-1):m*i,InvMahalanobis(Jx, MM_cov));
                [I_s2, J_s2, V_s2] = sparseFormat(m*i-(m-1):m*i,m*i+1:m*i+m,InvMahalanobis(-eye(m), MM_cov));
                I(tmp*i-(tmp-1):tmp*i) = horzcat(I_s1, I_s2); J(tmp*i-(tmp-1):tmp*i) = horzcat(J_s1, J_s2); V(tmp*i-(tmp-1):tmp*i) = horzcat(V_s1, V_s2);
                        
                X_diff = -X_next + X_pred;
                MM_vec(m*i-(m-1):m*i) = InvMahalanobis(X_diff, MM_cov);
            end
            MM_block = sparse(I, J, V, blk_height, blk_width);
            dummy = sparse([],[],[],blk_height,1);
            
            obj.opt.MM_vec = MM_vec;
            obj.opt.MM_block = horzcat(dummy, MM_block);

        end
        
        %% Wheel Speed Measurement Block and Residual
        function obj = CreateWSBlock(obj, states, wsf)
            m = size(states,1);
            n = size(states,2);
            
            blk_height = 3*n;

            if strcmp(obj.mode,'full')
                blk_width = 1 + m*n + 2*n*obj.prev_num;
            elseif strcmp(obj.mode,'partial')
                blk_width = 1 + m*n;
            end
            
            WS_vec = zeros(blk_height,1);

            cov = 5e-4*eye(3);
            tmp = 3+3*m;
            I = zeros(1,tmp*n); J = I; V = I;
            
            for i=1:n
                X = states(:,i);
                [Jwsf, Jx] = JacobianWS(X,wsf,obj.lr,obj.lane.xs.w(3,obj.dstart+i-1));
                [Iw, Jw, Vw] = sparseFormat(3*i-2:3*i,1,InvMahalanobis(Jwsf,cov));
                [Iv, Jv, Vv] = sparseFormat(3*i-2:3*i,1+m*i-(m-1):1+m*i,InvMahalanobis(Jx,cov));

                I(tmp*i-(tmp-1):tmp*i) = horzcat(Iw,Iv);
                J(tmp*i-(tmp-1):tmp*i) = horzcat(Jw,Jv);
                V(tmp*i-(tmp-1):tmp*i) = horzcat(Vw,Vv);

                V_b = VelMeas(X,wsf,obj.lr,obj.lane.xs.w(3,obj.dstart+i-1));
                % V_b_meas = [sqrt(obj.lane.vel.x(obj.dstart+i-1)^2 + obj.lane.vel.y(obj.dstart+i-1)^2);0;0];
                r = [obj.lane.eul.x(obj.dstart+i-1); obj.lane.eul.y(obj.dstart+i-1); obj.lane.eul.z(obj.dstart+i-1)];
                R = Euler_to_CTM(r);
                v = [obj.lane.vel.x(obj.dstart+i-1); obj.lane.vel.y(obj.dstart+i-1); obj.lane.vel.z(obj.dstart+i-1)];
                v_b = R * v;
                V_b_meas = [v_b(1); 0; 0];
                V_b_diff = -V_b_meas + V_b;

                WS_vec(3*i-2:3*i) = InvMahalanobis(V_b_diff,cov);
            end


            WS_block = sparse(I,J,V,blk_height,blk_width);

            obj.opt.WS_vec = WS_vec;
            obj.opt.WS_block = WS_block;
        end

        %% Measurement Model Block and Residual 
        function obj = CreateMEBlock(obj, states, latL, latR)
            
            if strcmp(obj.mode,'full') 
                [m,n] = size(states);
                x = 0:10:10*(obj.prev_num-1);
                blk_height = sum(obj.opt.assoc_idxsL(:,1:obj.prev_num),'all') + ...
                             sum(obj.opt.assoc_idxsR(:,1:obj.prev_num),'all');
                blk_width = m*n + 2*obj.prev_num*n;
                cnt = 1;
                ME_vec = zeros(blk_height,1);
                I = zeros(1,7*blk_height); J = I; V = I;
                for i=2:n
                    Xi = states(:,i);

                    for j=1:obj.prev_num
                        k_max_l = obj.opt.assoc_idxsL(i,j);
                        k_max_r = obj.opt.assoc_idxsR(i,j);

                        % Left Lane
                        for k=1:k_max_l
                            Xt = states(:,i-k);
                            lt_l = latL(i-k,j);               

                            [Jxt,Jxi,Jlt_l] = JacobianME2(Xt,Xi,lt_l,j,latL,i);
                            x_rel = Zpred_l(1); y_rel = Zpred_l(2);
                            
                            
                            
                            Lalp = latL(i,indic);
                            Lalpp1 = latL(i,indic+1);
                            

                            y = LatL(i,1:obj.prev_num);
                            cov = obj.output.cam.lstd_inter(obj.dstart+i-1,1:obj.prev_num);
                            y_meas = interp1(x,y,x_rel,'extrap');
                            cov_l = interp1(x,cov,x_rel,'extrap')^2;
                            
                            diff_l = -y_meas + y_rel;
                            ME_vec(cnt) = InvMahalanobis(diff_l,cov_l); 

                            [I1,J1,V1] = sparseFormat(cnt,m*(i-1)+1:m*(i-1)+2,InvMahalanobis(Jxi(1:2),cov_l));
                            [I2,J2,V2] = sparseFormat(cnt,m*(i-1)+9,InvMahalanobis(Jxi(3),cov_l));
                            [I3,J3,V3] = sparseFormat(cnt,m*(i-k-1)+1:m*(i-k-1)+2,InvMahalanobis(Jxt(1:2),cov_l));
                            [I4,J4,V4] = sparseFormat(cnt,m*(i-k-1)+9,InvMahalanobis(Jxt(3),cov_l));
                            [I5,J5,V5] = sparseFormat(cnt,m*n+(i-k-1)*obj.prev_num+j,InvMahalanobis(Jlt_l,cov_l));

                            I(7*cnt-6:7*cnt) = horzcat(I1,I2,I3,I4,I5);
                            J(7*cnt-6:7*cnt) = horzcat(J1,J2,J3,J4,J5);
                            V(7*cnt-6:7*cnt) = horzcat(V1,V2,V3,V4,V5);

                            cnt = cnt + 1;
                        end
                        % Right Lane
                        for k=1:k_max_r
                            Xt = states(:,i-k);
                            lt_r = latR(i-k,j);
                            [Jxt,Jxi,Jlt_r,Zpred_r] = JacobianME2(Xt,Xi,lt_r,j);

                            x_rel = Zpred_r(1); y_rel = Zpred_r(2);
                            y = obj.output.cam.r_inter(obj.dstart+i-1,1:obj.prev_num);
                            cov = obj.output.cam.rstd_inter(obj.dstart+i-1,1:obj.prev_num);
                            y_meas = interp1(x,y,x_rel);
                            cov_r = interp1(x,cov,x_rel)^2;

                            diff_r = -y_meas + y_rel;
                            ME_vec(cnt) = InvMahalanobis(diff_r,cov_r); 
                            [I1,J1,V1] = sparseFormat(cnt,m*(i-1)+1:m*(i-1)+2,InvMahalanobis(Jxi(1:2),cov_r));
                            [I2,J2,V2] = sparseFormat(cnt,m*(i-1)+9,InvMahalanobis(Jxi(3),cov_r));
                            [I3,J3,V3] = sparseFormat(cnt,m*(i-k-1)+1:m*(i-k-1)+2,InvMahalanobis(Jxt(1:2),cov_r));
                            [I4,J4,V4] = sparseFormat(cnt,m*(i-k-1)+9,InvMahalanobis(Jxt(3),cov_r));
                            [I5,J5,V5] = sparseFormat(cnt,m*n+(n+i-k-1)*obj.prev_num+j,InvMahalanobis(Jlt_r,cov_r));

                            I(7*cnt-6:7*cnt) = horzcat(I1,I2,I3,I4,I5);
                            J(7*cnt-6:7*cnt) = horzcat(J1,J2,J3,J4,J5);
                            V(7*cnt-6:7*cnt) = horzcat(V1,V2,V3,V4,V5);

                            cnt = cnt + 1;
                        end
                    end
                end
                dummy = sparse([],[],[],blk_height,1);
                ME_block = horzcat(dummy, sparse(I,J,V,blk_height,blk_width));              
                obj.opt.ME_block = ME_block;
                obj.opt.ME_vec = ME_vec;

            elseif strcmp(obj.mode,'partial')
                obj.opt.ME_vec = [];
                obj.opt.ME_block = [];
            end
        end

        %% Merge Block and Residuals
        function obj = Merge(obj)
        % Merge multi-modal residuals and jacobian block matrices
        % Passed to NLS solver for optimization

            obj.opt.res = vertcat(obj.opt.Pr_vec, ...
                                  obj.opt.GNSS_vec, ...
                                  obj.opt.MM_vec, ...
                                  obj.opt.WS_vec, ...
                                  obj.opt.ME_vec);

            obj.opt.jac = vertcat(obj.opt.Pr_block, ...
                                  obj.opt.GNSS_block, ...
                                  obj.opt.MM_block, ...
                                  obj.opt.WS_block, ...
                                  obj.opt.ME_block);
        end

        %% Compute 2D Lateral and Longitudinal Error
        function obj = computeError(obj)
            disp('-Computing Optimization Error-')
            gtx = obj.lane.posx(obj.dstart:obj.dend);
            gty = obj.lane.posy(obj.dstart:obj.dend);
            gtpsi = obj.lane.eul.z(obj.dstart:obj.dend);
            
            x = obj.opt.states(1,:);
            y = obj.opt.states(2,:);
            
            n = size(x,2); % number of states
            obj.opt.error = struct();
            obj.opt.error.lon = zeros(1,n);
            obj.opt.error.lat = zeros(1,n);
            obj.opt.error.lon_mse = 0; obj.opt.error.lat_mse = 0; obj.opt.error.mse = 0;
            obj.opt.error.e = zeros(1,n);
            
            for i=1:n
                delx = x(i) - gtx(i); dely = y(i) - gty(i);
                obj.opt.error.e(i) = sqrt(delx^2 + dely^2);
                obj.opt.error.mse = obj.opt.error.mse + delx^2 + dely^2;
                obj.opt.error.lon(i) = delx * cos(gtpsi(i)) + dely * sin(gtpsi(i));
                obj.opt.error.lat(i) = -delx * sin(gtpsi(i)) + dely * cos(gtpsi(i));
                obj.opt.error.lon_mse = obj.opt.error.lon_mse + obj.opt.error.lon(i)^2;
                obj.opt.error.lat_mse = obj.opt.error.lat_mse + obj.opt.error.lat(i)^2;
            end
            
            obj.opt.error.rmse = sqrt(obj.opt.error.mse/n);
            obj.opt.error.lon_rmse = sqrt(obj.opt.error.lon_mse/n);
            obj.opt.error.lat_rmse = sqrt(obj.opt.error.lat_mse/n);
            obj.opt.error.maxe = max(abs(obj.opt.error.e));
            obj.opt.error.maxlon = max(abs(obj.opt.error.lon));
            obj.opt.error.maxlat = max(abs(obj.opt.error.lat));
            disp('==================================================')
            disp(['Optimization RMSE: ', num2str(obj.opt.error.rmse), 'm'])
            disp(['Optimization Longitudinal RMSE: ', num2str(obj.opt.error.lon_rmse), 'm'])
            disp(['Optimization Lateral RMSE: ', num2str(obj.opt.error.lat_rmse), 'm'])
            disp(['Maximum Longitudinal Error: ', num2str(obj.opt.error.maxlon), 'm'])
            disp(['Maximum Lateral Error: ', num2str(obj.opt.error.maxlat), 'm'])
            disp(['Maximum Error: ', num2str(obj.opt.error.maxe), 'm'])
            disp('==================================================')
        end
        
        %% Plot Results
        function obj = plotRes(obj)
            disp('-Plotting Optimization Results-')
            
            % 2D Localization and Mapping Result
            mapbox_opts.accesstoken = 'sk.eyJ1IjoiamluaHdhbjk4IiwiYSI6ImNreTU3dWhjMzBqZDkydm9ud3MyeWNiaWwifQ.ViL0qz-v5MF9Ggxym4AwIQ';
            plot_opts = struct();
            plot_opts.maptype = 'satellite';
            plot_opts.scale = 1;
            plot_opts.resize = 1;

            figure(1);

            ref_lat = mean(obj.lane.lat); ref_lon = mean(obj.lane.lon);
            ref = [ref_lat; ref_lon];
            x = obj.opt.states(1,:); y = obj.opt.states(2,:);
            gnss_x = obj.gnss.meas(1,:); gnss_y = obj.gnss.meas(2,:);
            
            ins_x = obj.ins.state_saver(1,:); ins_y = obj.ins.state_saver(2,:);
            [ins_lon, ins_lat] = lin_to_geo(ins_x, ins_y, ref);

            % Plot 0m previewed optimized lane points
            [opt_veh_lon, opt_veh_lat] = lin_to_geo(x, y, ref);
            [gnss_lon, gnss_lat] = lin_to_geo(gnss_x, gnss_y, ref);

            p_gnss = plot(gnss_x, gnss_y,'gx'); hold on; axis equal; grid on;

%             p_gnss = plot(gnss_lon, gnss_lat,'Color',[1 165/255 0],'Marker','.','LineStyle','none'); hold on; axis equal; grid on;
            p_opt2 = plot(x, y,'b.'); 
            p_opt1 = plot(x(obj.gnss.idxs-obj.dstart+1), y(obj.gnss.idxs-obj.dstart+1),'r.');
            
%             p_ins = plot(ins_lon, ins_lat,'m.');


            [m,n] = size(obj.opt.states);
            if strcmp(obj.mode,'full') 
                                
                for i=1:obj.prev_num
                    opt_lml = reshape(obj.opt.lml(:,i),2,[]);
                    opt_lmr = reshape(obj.opt.lmr(:,i),2,[]);
                                       
                    [opt_lml_lon, opt_lml_lat] = lin_to_geo(opt_lml(1,:), opt_lml(2,:), ref);
                    [opt_lmr_lon, opt_lmr_lat] = lin_to_geo(opt_lmr(1,:), opt_lmr(2,:), ref);
    
                    p_left = plot(opt_lml(1,:), opt_lml(2,:));
                    p_right = plot(opt_lmr(1,:), opt_lmr(2,:));

%                     ins_lml = reshape(obj.ins.lml_saver(:,i),2,[]);
%                     ins_lmr = reshape(obj.ins.lmr_saver(:,i),2,[]);
%                     [ins_lml_lon, ins_lml_lat] = lin_to_geo(ins_lml(1,:), ins_lml(2,:), ref);
%                     [ins_lmr_lon, ins_lmr_lat] = lin_to_geo(ins_lmr(1,:), ins_lmr(2,:), ref);
% 
%                     plot(ins_lml_lon, ins_lml_lat,'m--');
%                     plot(ins_lmr_lon, ins_lmr_lat,'c--');
                end
            end

            p_gt = plot(obj.lane.posx(obj.dstart:obj.dend), obj.lane.posy(obj.dstart:obj.dend),'k--');

%             plot_mapbox('maptype',plot_opts.maptype,'scale',plot_opts.scale,'resize',plot_opts.resize,'apikey',mapbox_opts.accesstoken)

            xlabel('Longitude(Deg)'); ylabel('Latitude(Deg)'); 
            
            
            if strcmp(obj.mode,'full')
                title('GNSS + INS + WSS + Lane Detection Sensor Fusion via MAP');
                legend([p_opt1 p_opt2 p_gt p_gnss p_left p_right],...
                       'Optimized Trajectory(w GNSS)','Optimized Trajectory(w/o GNSS)',...
                       'Ground Truth','GNSS Measurement(Ublox)',...
                       'Optimized Left Lane', 'Optimized Right Lane');
            
            elseif strcmp(obj.mode,'partial')
                title('GNSS + INS + WSS Sensor Fusion via MAP');
                legend([p_opt1, p_opt2, p_gt, p_gnss],...
                       'Optimized Trajectory(w GNSS)','Optimized Trajectory(w/o GNSS)',...
                       'Ground Truth','GNSS Measurement(Ublox)')
            end

            t = obj.lane.t(obj.dstart:obj.dend) - obj.lane.t(obj.dstart);
            %% Vehicle Heading Angle
            figure(2);
            subplot(3,1,1)
            phi_diff = 180/pi * obj.opt.states(7,:) - 180/pi * obj.lane.eul.x(obj.dstart:obj.dend)';
            p_diff2 = plot(t,phi_diff,'b.'); hold on; grid on; axis tight;
            p_diff1 = plot(t(obj.gnss.idxs-obj.dstart+1),phi_diff(obj.gnss.idxs-obj.dstart+1),'r.'); 

            xlabel('Time(s)'); ylabel('Orientation(Deg)'); title('Roll Angle Difference w.r.t Ground Truth');
            legend([p_diff1, p_diff2],'Optimized(w GNSS)', 'Optimized(w/o GNSS)')

            subplot(3,1,2)
            theta_diff = 180/pi * obj.opt.states(8,:) - 180/pi * obj.lane.eul.y(obj.dstart:obj.dend)';
            p_diff2 = plot(t,theta_diff,'b.'); hold on; grid on; axis tight;
            p_diff1 = plot(t(obj.gnss.idxs-obj.dstart+1),theta_diff(obj.gnss.idxs-obj.dstart+1),'r.'); 

            xlabel('Time(s)'); ylabel('Orientation(Deg)'); title('Pitch Angle Difference w.r.t Ground Truth');
            legend([p_diff1, p_diff2],'Optimized(w GNSS)', 'Optimized(w/o GNSS)')

            subplot(3,1,3)
            psi_diff = 180/pi * obj.opt.states(9,:) - 180/pi * obj.lane.eul.z(obj.dstart:obj.dend)';
            p_diff2 = plot(t,psi_diff,'b.'); hold on; grid on; axis tight;
            p_diff1 = plot(t(obj.gnss.idxs-obj.dstart+1),psi_diff(obj.gnss.idxs-obj.dstart+1),'r.'); 

            xlabel('Time(s)'); ylabel('Orientation(Deg)'); title('Yaw Angle Difference w.r.t Ground Truth');
            legend([p_diff1, p_diff2],'Optimized(w GNSS)', 'Optimized(w/o GNSS)')

            %% Vehicle Velocity Comparison
            figure(3);
            subplot(3,1,1);
            vx = obj.opt.states(4,:);
            p_opt2 = plot(t,vx,'b.'); hold on; grid on; axis tight;
            p_opt1 = plot(t(obj.gnss.idxs-obj.dstart+1),vx(obj.gnss.idxs-obj.dstart+1),'r.'); 

            p_gt = plot(t,obj.lane.vel.x(obj.dstart:obj.dend));
            xlabel('Time(s)'); ylabel('Vx (m/s)'); title('Vx Comparison');
            legend([p_opt1, p_opt2, p_gt],'Optimized(w GNSS)','Optimized(w/o GNSS)','Ground Truth')

            subplot(3,1,2);
            vy = obj.opt.states(5,:);
            p_opt2 = plot(t,vy,'b.'); hold on; grid on; axis tight;
            p_opt1 = plot(t(obj.gnss.idxs-obj.dstart+1),vy(obj.gnss.idxs-obj.dstart+1),'r.'); 

            p_gt = plot(t,obj.lane.vel.y(obj.dstart:obj.dend));
            xlabel('Time(s)'); ylabel('Vy (m/s)'); title('Vy Comparison');
            legend([p_opt1, p_opt2, p_gt],'Optimized(w GNSS)','Optimized(w/o GNSS)','Ground Truth')

            subplot(3,1,3);
            vz = obj.opt.states(6,:);
            p_opt2 = plot(t,vz,'b.'); hold on; grid on; axis tight;
            p_opt1 = plot(t(obj.gnss.idxs-obj.dstart+1),vz(obj.gnss.idxs-obj.dstart+1),'r.'); 

            p_gt = plot(t,obj.lane.vel.z(obj.dstart:obj.dend));
            xlabel('Time(s)'); ylabel('Vz (m/s)'); title('Vz Comparison');
            legend([p_opt1, p_opt2, p_gt],'Optimized(w GNSS)','Optimized(w/o GNSS)','Ground Truth')

            %% Euler Angle Comparison
            figure(4);
            subplot(3,1,1);
            eul_x = 180/pi* obj.opt.states(7,:);
            p_opt2 = plot(t,eul_x,'b.'); hold on; grid on; axis tight;
            p_opt1 = plot(t(obj.gnss.idxs-obj.dstart+1),eul_x(obj.gnss.idxs-obj.dstart+1),'r.'); 

            p_gt = plot(t,180/pi* obj.lane.eul.x(obj.dstart:obj.dend));
            xlabel('Time(s)'); ylabel('Euler X (Deg)'); title('Euler X Comparison');
            legend([p_opt1, p_opt2, p_gt],'Optimized(w GNSS)','Optimized(w/o GNSS)','Ground Truth')

            subplot(3,1,2);
            eul_y = 180/pi* obj.opt.states(8,:);
            p_opt2 = plot(t,eul_y,'b.'); hold on; grid on; axis tight;
            p_opt1 = plot(t(obj.gnss.idxs-obj.dstart+1),eul_y(obj.gnss.idxs-obj.dstart+1),'r.'); 

            p_gt = plot(t,180/pi* obj.lane.eul.y(obj.dstart:obj.dend));
            xlabel('Time(s)'); ylabel('Euler Y (Deg)'); title('Euler Y Comparison');
            legend([p_opt1, p_opt2, p_gt],'Optimized(w GNSS)','Optimized(w/o GNSS)','Ground Truth')

            subplot(3,1,3);
            eul_z = 180/pi* obj.opt.states(9,:);
            p_opt2 = plot(t,eul_z,'b.'); hold on; grid on; axis tight;
            p_opt1 = plot(t(obj.gnss.idxs-obj.dstart+1),eul_z(obj.gnss.idxs-obj.dstart+1),'r.'); 

            p_gt = plot(t,180/pi* obj.lane.eul.z(obj.dstart:obj.dend));
            xlabel('Time(s)'); ylabel('Euler Z (Deg)'); title('Euler Z Comparison');
            legend([p_opt1, p_opt2, p_gt],'Optimized(w GNSS)','Optimized(w/o GNSS)','Ground Truth')

            %% Acceleration Bias
            figure(5);

            abx = obj.opt.states(10,:);
            aby = obj.opt.states(11,:);
            abz = obj.opt.states(12,:);
            p_abx = plot(t,abx,'r.'); hold on; grid on; axis tight;
            p_aby = plot(t,aby,'g.');
            p_abz = plot(t,abz,'b.');
            xlabel('Time(s)'); ylabel('Bias(m/s^2)'); title('Acceleration Bias');
            legend([p_abx, p_aby, p_abz],'Abx','Aby','Abz')

            %% Acceleration Bias
            figure(6);

            wbx = 180/pi * obj.opt.states(13,:);
            wby = 180/pi * obj.opt.states(14,:);
            wbz = 180/pi * obj.opt.states(15,:);
            p_wbx = plot(t,wbx,'r.'); hold on; grid on; axis tight;
            p_wby = plot(t,wby,'g.');
            p_wbz = plot(t,wbz,'b.');
            xlabel('Time(s)'); ylabel('Bias(deg/s)'); title('Angular Velocity Bias');
            legend([p_wbx, p_wby, p_wbz],'Wbx','Wby','Wbz')

            %% Compute RMSE

            elon = obj.opt.error.lon; elat = obj.opt.error.lat;
                    
            figure(7);
            subplot(2,1,1);
            plon1 = plot(t,elon,'b.'); hold on; grid on;
            plon2 = plot(t(obj.gnss.idxs-obj.dstart+1),elon(obj.gnss.idxs-obj.dstart+1),'r.');
            xlabel('Time(s)'); ylabel('Error(m)'); title('Longitudinal Error');
            legend([plon1, plon2],'Error(w GNSS)','Error(w/o GNSS)')

            subplot(2,1,2);
            plat1 = plot(t,elat,'b.'); hold on; grid on;
            plat2 = plot(t(obj.gnss.idxs-obj.dstart+1),elat(obj.gnss.idxs-obj.dstart+1),'r.');
            xlabel('Time(s)'); ylabel('Error(m)'); title('Lateral Error');
            legend([plat1, plat2],'Error(w GNSS)','Error(w/o GNSS)')

            figure(8);
            e = obj.opt.error.e;
            e1 = plot(t,e,'b.'); hold on; grid on;
            e2 = plot(t(obj.gnss.idxs-obj.dstart+1), e(obj.gnss.idxs-obj.dstart+1),'r.');
            xlabel('Time(s)'); ylabel('Error(m)'); title('2D Error')
            legend([e2, e1],'Error(w GNSS)','Error(w/o GNSS)')

        end
        
        %% Re-order lane points
        function obj = reorderLP(obj)
            
            % Need to modify accessing covariance part
            disp('-Lane Point Re-ordering for Parametrization-')
            [m,n] = size(obj.ins.state_saver);
            
            lml_pc = reshape(reshapeLP(obj.opt.lml),2,[]);
            lmr_pc = reshape(reshapeLP(obj.opt.lmr),2,[]);
            
            lml_cov = extractCov(obj.opt.cov(2+m*n:1+m*n+2*n*obj.prev_num,2+m*n:1+m*n+2*n*obj.prev_num));
            lmr_cov = extractCov(obj.opt.cov(2+m*n+2*n*obj.prev_num:end,2+m*n+2*n*obj.prev_num:end));           

            [obj.opt.reordered_lml_pc, obj.opt.reordered_lmr_pc, ...
             obj.opt.reordered_lml_cov, obj.opt.reordered_lmr_cov] = sortLP(lml_pc, lmr_pc, lml_cov, lmr_cov);
        end

        %% Compute Measurement Jacobian
        function [Jxt,Jxo,Jlt,JLalp,JLalpp1,resi] = JacobianME2(obj,X_t,X_o,l_t,idx,lat,i)
            x_t = X_t(1); y_t = X_t(2); psi_t = X_t(9);
            x_o = X_o(1); y_o = X_o(2); psi_o = X_o(9);
            L_tj  = [x_t + 10*(idx-1)*cos(psi_t) - l_t*sin(psi_t);
                     y_t + 10*(idx-1)*sin(psi_t) + l_t*cos(psi_t)];
            delX = L_tj - [x_o; y_o];
           
            Xb = [cos(psi_o) sin(psi_o);-sin(psi_o) cos(psi_o)] * delX;
            xb = Xb(1); yb = Xb(2);

            indic = floor(x_b/10)+1;  % Indicator for preview number lower bound     
            if indic < 1
                indic = 1;
            elseif indic >= obj.prev_num
                indic = obj.prev_num-1;
            end
            Lalp = lat(i,indic);
            Lalpp1 = lat(i)
            ymeas = (Lalpp1 - Lalp)/10 * (xb - 10*(indic-1)) + Lalp;
        
            resi = yb - ymeas;
            Jxt = zeros(1,3);
            Jxt(1) = -sin(psi_o); Jxt(2) = cos(psi_o); 
            Jxt(3) = cos(psi_o)*(cos(psi_t)*(10*idx - 10) - l_t*sin(psi_t)) + sin(psi_o)*(sin(psi_t)*(10*idx - 10) + l_t*cos(psi_t));
            
            Jxo = zeros(1,3);
            Jxo(1) = sin(psi_o); Jxo(2) = -cos(psi_o);
            Jxo(3) = cos(psi_o)*(x_o - x_t - cos(psi_t)*(10*idx - 10) + l_t*sin(psi_t)) - sin(psi_o)*(y_t - y_o + sin(psi_t)*(10*idx - 10) + l_t*cos(psi_t));
            
            Jlt = cos(psi_o - psi_t);
        
            JLalp = (sin(psi_o)*(y_t - y_o + sin(psi_t)*(10*idx - 10) + l_t*cos(psi_t)))/10 - indic - (cos(psi_o)*(x_o - x_t - cos(psi_t)*(10*idx - 10) + l_t*sin(psi_t)))/10;
            
            JLalpp1 = indic - (sin(psi_o)*(y_t - y_o + sin(psi_t)*(10*idx - 10) + l_t*cos(psi_t)))/10 + (cos(psi_o)*(x_o - x_t - cos(psi_t)*(10*idx - 10) + l_t*sin(psi_t)))/10 - 1;
        
        end
    end
end

%% ======================================= Other Functions =======================================
%% State Prediction
function X_next = StatePred(X_curr, u, dt)
    rx = X_curr(7); ry = X_curr(8); rz = X_curr(9);
    a = u(1:3); w = u(4:6);
    
    % Noise std
    na_std = 6.9e-4;
    nw_std = 1.49e-4;
    nad_cov = 1/dt * diag(ones(1,3) * na_std^2);
    nwd_cov = 1/dt * diag(ones(1,3) * nw_std^2);
    nad = mvnrnd(zeros(3,1), nad_cov)';
    nwd = mvnrnd(zeros(3,1), nwd_cov)';
    
    % Bias std
    ab_std = 0.9e-5;
    wb_std = 2*1e-06;
    abd_cov = 1/dt * diag(ones(1,3) * ab_std^2);
    wbd_cov = 1/dt * diag(ones(1,3) * wb_std^2);
    abd = mvnrnd(zeros(3,1), abd_cov)';
    wbd = mvnrnd(zeros(3,1), wbd_cov)';
    
    % Euler angle --> Rotation Matrix
    CTM = Euler_to_CTM([rx; ry; rz]);
    CTMbn = CTM';
    % Attitude Update
    wb = X_curr(13:15);
    w_n = w - wb - nwd;
    d_euler = w_n * dt;
    A = Exp_map(d_euler);
    
    CTMbn_n = CTMbn * A; 
    CTMnb_n = CTMbn_n';
    r_n = CTM_to_Euler(CTMnb_n);
    
    % Velocity Update
    v = X_curr(4:6); ab = X_curr(10:12);
    g_ = [0; 0; -9.81];
    
    v_n = v + g_ * dt + CTMbn * (a - ab - nad) * dt;
    
    % Position Update
    p = X_curr(1:3);
    p_n = p + v * dt + 1/2 * g_ * dt^2 + 1/2 * CTMbn * (a - ab - nad) * dt^2;
    
    % Bias Update
    ab_n = ab + abd * dt;
    wb_n = wb + wbd * dt;
    
    X_next = vertcat(p_n, v_n, r_n, ab_n, wb_n);
end

%% Lane Point Prediction
% function L_next = LanePred(X_next, Lij, Lijp1, idx)
%     
%     L_next = zeros(2,1);
%     
%     x = X_next(1); y = X_next(2); psi = X_next(9);
%     lijp1x = Lijp1(1); lijp1y = Lijp1(2);
%     lijx = Lij(1); lijy = Lij(2);
%     L_next(1) = x + cos(psi)*(10*idx - 10) - (sin(psi)*(lijy - y - sin(psi)*(10*idx - 10) + ((lijy - lijp1y)*(x - lijx + cos(psi)*(10*idx - 10)))/(lijx - lijp1x)))/(cos(psi) + (sin(psi)*(lijy - lijp1y))/(lijx - lijp1x));
%     L_next(2) = y + sin(psi)*(10*idx - 10) + (cos(psi)*(lijy - y - sin(psi)*(10*idx - 10) + ((lijy - lijp1y)*(x - lijx + cos(psi)*(10*idx - 10)))/(lijx - lijp1x)))/(cos(psi) + (sin(psi)*(lijy - lijp1y))/(lijx - lijp1x));
% end

%% Body Frame velocity computation
function V_b = VelMeas(X,wsf,ls,wz)
    % Convert world frame velocity values to vehicle body frame velocity
    CTMnb = Euler_to_CTM(X(7:9));
%     - [0;ls*(wz-X(15));0]
    V_b = 1/wsf * (CTMnb * X(4:6) - [0;ls*(wz-X(15));0]); 
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

%% Latitude/Longitude to TM Coord. Converter
function [PosLon,PosLat]=lin_to_geo(localX, localY, PosRef)

% PosLat : Latitude in Degrees
% PosLon : Longitude in Degrees
% PosRef : Origin of Local coordinate
% lpos : [ localX, localY ] in meter, East-North Coordinate

% Convert Geographic coordinate into Linear coordinate with WGS84
% Ellipsoid model constants (actual values here are for WGS84

R0 = 6378137.0;
E=1/298.257223563;

Rn = R0*(1-E^2)/((1-E^2 * (sind(PosRef(1))^2))^(3/2));
Re = R0/((1-E^2 *(sind(PosRef(1))^2))^(1/2));
Ra = Re*Rn / sqrt( Re^2*sind(PosRef(1))^2 + Rn^2*cosd(PosRef(1))^2 );

deltaLon = localX * 180/pi /Ra / cosd(PosRef(1));
deltaLat = localY * 180/pi /Rn;

PosLon = deltaLon + PosRef(2);
PosLat = deltaLat + PosRef(1);

end

%% Get Lane Point from lateral info
function LP = getLP(X,lat,idx)
    LP = zeros(2,1);
    x = X(1); y = X(2); psi = X(9);
    LP(1) = x + 10*(idx-1)*cos(psi) - lat*sin(psi);
    LP(2) = y + 10*(idx-1)*sin(psi) + lat*cos(psi);
end

%% Reshape Lane Point for creating intial guess
function res = reshapeLP(LP)
    num_states = size(LP,1)/2;
    prev_num = size(LP,2);
    
    res = zeros(size(LP,1)*size(LP,2),1);
    
    for i=1:num_states
        tmp = (2*prev_num)*i-(2*prev_num);
        for j=1:prev_num
            LP_sample = LP(2*i-1:2*i,j);
            res(tmp + 2*j-1:tmp + 2*j) = LP_sample;
        end
    end
end

%% Extract covariance from information matrix
function colCov = extractCov(cov_mat)
    n = size(cov_mat,1)/2;
    colCov = zeros(4,n);
    for i=1:n
        cov = cov_mat(2*i-1:2*i,2*i-1:2*i);
        colCov(:,i) = reshape(cov,4,1);
    end
end

%% Sort LP using nearest Neighbor
function [pc_re_l, pc_re_r, cov_re_l, cov_re_r] = sortLP(pc_l, pc_r, covs_l, covs_r)
    search_idx = 1; 
    n = size(pc_l,2);
    
    pc_lcpyd = vertcat(pc_l(:,2:end), 2:1:n); % Add absolute indices
    pc_rcpyd = vertcat(pc_r(:,2:end), 2:1:n);
    cov_lcpyd = covs_l(:,2:end);
    cov_rcpyd = covs_r(:,2:end);
    pc_l = vertcat(pc_l, 1:n);
    pc_r = vertcat(pc_r, 1:n);
    
    pc_re_l = pc_l(:,1);
    pc_re_r = pc_r(:,1);
    cov_re_l = covs_l(:,1);
    cov_re_r = covs_r(:,1);

    % Phase 1: Sort Lane Points by iteratively finding the nearest point
    disp('Phase1: Nearest Neighbor Search')
    while search_idx ~= n
        px_l = pc_l(1,search_idx); py_l = pc_l(2,search_idx);
        px_r = pc_r(1,search_idx); py_r = pc_r(2,search_idx);

        d = (pc_lcpyd(1,:) - px_l).^2 + (pc_lcpyd(2,:) - py_l).^2 + ...
            (pc_rcpyd(1,:) - px_r).^2 + (pc_rcpyd(2,:) - py_r).^2;
        
        [~,idx] = min(d);
        search_idx = pc_lcpyd(3,idx);
        pc_lcpyd(:,idx) = []; pc_rcpyd(:,idx) = [];
        cov_lcpyd(:,idx) = []; cov_rcpyd(:,idx) = [];

        pc_re_l = [pc_re_l pc_l(:,search_idx)];
        pc_re_r = [pc_re_r pc_r(:,search_idx)];
        cov_re_l = [cov_re_l covs_l(:,search_idx)];
        cov_re_r = [cov_re_r covs_r(:,search_idx)];
    end

    left_pc_l = pc_lcpyd; left_pc_r = pc_rcpyd;
    left_cov_l = cov_lcpyd; left_cov_r = cov_rcpyd;
    
    % Phase 2: Fill in left points
    disp('Phase2: Fill in left out points')
    n = size(left_pc_l,2);
    for i=1:n
        px_l = left_pc_l(1,i); py_l = left_pc_l(2,i);
        px_r = left_pc_r(1,i); py_r = left_pc_r(2,i);

        d = (pc_re_l(1,:) - px_l).^2 + (pc_re_l(2,:) - py_l).^2 + ...
            (pc_re_r(1,:) - px_r).^2 + (pc_re_l(2,:) - py_r).^2;
        [~,idx] = min(d);

        if d(idx-1) > d(idx+1)
            pc_re_l = horzcat(pc_re_l(:,1:idx), left_pc_l(:,i), pc_re_l(:,idx+1:end));
            pc_re_r = horzcat(pc_re_r(:,1:idx), left_pc_r(:,i), pc_re_r(:,idx+1:end));
            cov_re_l = horzcat(cov_re_l(:,1:idx), left_cov_l(:,i), cov_re_l(:,idx+1:end));
            cov_re_r = horzcat(cov_re_r(:,1:idx), left_cov_r(:,i), cov_re_r(:,idx+1:end));
        else
            pc_re_l = horzcat(pc_re_l(:,1:idx-1), left_pc_l(:,i), pc_re_l(:,idx:end));
            pc_re_r = horzcat(pc_re_r(:,1:idx-1), left_pc_r(:,i), pc_re_r(:,idx:end));
            cov_re_l = horzcat(cov_re_l(:,1:idx), left_cov_l(:,i), cov_re_l(:,idx+1:end));
            cov_re_r = horzcat(cov_re_r(:,1:idx), left_cov_r(:,i), cov_re_r(:,idx+1:end));
        end
    end

    % Phase 3: Perform Circular Fitting and finish re-ordering
    disp('Phase 3: Perform Circular Fitting to accurately order points')
    
    fit_length = 20;
    n = size(pc_re_l,2);
    for i=1:n-fit_length+1
        disp(num2str('Iteration: ',num2str(i)))
        [res,~] = Circle2DFitV2(pc_re_l,pc_re_r,cov_re_l,cov_re_r,0.9,[i i+fit_length-1],false);
        if res.th(1) < res.th(end)
            [~,sorted_idxs] = sort(res.th);
        else
            [~,sorted_idxs] = sort(res.th,'descend');
        end
        disp(sorted_idxs)
        
        pcs_l = pc_re_l(:,i:i+fit_length-1);
        pcs_r = pc_re_r(:,i:i+fit_length-1);
        covs_l = cov_re_l(:,i:i+fit_length-1);
        covs_r = cov_re_r(:,i:i+fit_length-1);
        
        pc_re_l(:,i:i+fit_length-1) = pcs_l(:,sorted_idxs);
        pc_re_r(:,i:i+fit_length-1) = pcs_r(:,sorted_idxs);
        cov_re_l(:,i:i+fit_length-1) = covs_l(:,sorted_idxs);
        cov_re_r(:,i:i+fit_length-1) = covs_r(:,sorted_idxs);
    end
end
