classdef Optimizer_LinV4 < handle
% Optimizer module for Batch Optimization
% Lane is represented as arcs
% Implemented by JinHwan Jeon, 2022
    properties
        lane 
        output
        prev_num % Preview distance index 1~10
        dstart % Dataset start index
        dend % Dataset end index
        x0 % Optimization solver intial value
        lr = 1.44; % Distance from CoG(RT location) to rear axle
        % ls = 0.8; % Distance from CoG to IMU, Camera (Sensor Lever Arm)
        ins = struct(); % INS Mechanization Saver
        gnss = struct(); % GNSS Measurement Saver
        opt = struct(); % Optimzation Results Saver        
        plot_flag % Flag for automatic result plot
    end
    methods
        %% Constructor
        function obj = Optimizer_LinV4(lane, output, prev_num, dstart, dend, plot_flag)
            obj.lane = lane;
            obj.output = output;
            obj.prev_num = prev_num;
            obj.dstart = dstart;
            obj.dend = dend;
            
            obj.plot_flag = plot_flag;

            obj.ins.dstart = dstart;
            obj.ins.dend = dend;
            obj.ins.iter = dstart;
            obj.ins.state = [lane.posx(dstart); lane.posy(dstart); lane.posz(dstart); 
                             lane.vel.x(dstart); lane.vel.y(dstart); lane.vel.z(dstart);
                             lane.eul.x(dstart); lane.eul.y(dstart); lane.eul.z(dstart);
                             0; 0; 0;
                             0; 0; 0];
            
            init_thetaL = atan2(lane.ly(dstart,2) - lane.ly(dstart,1), lane.lx(dstart,2) - lane.lx(dstart,1));
            init_kappaL = getCircle([lane.lx(dstart,1:3);lane.ly(dstart,1:3)]);

            init_thetaR = atan2(lane.ly(dstart,2) - lane.ly(dstart,1), lane.lx(dstart,2) - lane.lx(dstart,1));
            init_kappaR = getCircle([lane.lx(dstart,1:3);lane.ly(dstart,1:3)]);
            obj.ins.paramsL = [lane.lx(dstart,1); lane.ly(dstart,1); init_thetaL; init_kappaL];
            obj.ins.paramsR = [lane.rx(dstart,1); lane.ry(dstart,1); init_thetaR; init_kappaR];
            

            obj.ins.state_saver = zeros(size(obj.ins.state,1),dend-dstart+1);
            obj.ins.state_saver(:,1) = obj.ins.state;
           
            obj.opt.lb = []; obj.opt.ub = [];
            obj.opt.options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',false,'Display','iter-detailed','Algorithm','trust-region-reflective','UseParallel',true);
            obj.opt.states = [];            
            
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
            
            obj.x0 = vertcat(1,reshape(obj.ins.state_saver,[],1),...
                             obj.ins.paramsL,obj.ins.paramsR);                        
        end

        %% Pre-process GNSS measurements
        function obj = GNSS(obj)
            disp('-GNSS Pre-processing-')
            obj.gnss.lb = find(obj.output.ubloxgps.filtered_idxs == obj.dstart);
            obj.gnss.ub = find(obj.output.ubloxgps.filtered_idxs == obj.dend);
            
            if isempty(obj.gnss.lb) || isempty(obj.gnss.ub)
                error('No GNSS measurement at start or end points, select other data index for start or end points')
            end
            % scfull, sc1, sc2(high curvature)
%             if obj.output.ubloxgps.filtered_idxs(obj.gnss.lb) ~= obj.output.ubloxgps.sc2.trimmed_idxs(1) || ...
%                obj.output.ubloxgps.filtered_idxs(obj.gnss.ub) ~= obj.output.ubloxgps.sc2.trimmed_idxs(end)
%                 error('Trimmed GNSS indices need to be updated, default is 17301~22000')
%                 % To run optimization for different intervals, need to
%                 % correct "ubloxgps_proc.m" in Functions folder and save
%                 % updated results to "output.mat"
%             end

            obj.gnss.cov = obj.output.ubloxgps.filtered_cov(:,obj.gnss.lb:obj.gnss.ub);
            
            obj.gnss.idxs = obj.output.ubloxgps.filtered_idxs(obj.gnss.lb:obj.gnss.ub);
            % obj.gnss.iter_idxs = obj.gnss.idxs(obj.gnss.idxs <= obj.dstart);
            
            obj.gnss.meas = [obj.output.ubloxgps.filtered_x(obj.gnss.lb:obj.gnss.ub); 
                             obj.output.ubloxgps.filtered_y(obj.gnss.lb:obj.gnss.ub);
                             obj.output.ubloxgps.filtered_z(obj.gnss.lb:obj.gnss.ub)];   
            
%             obj.gnss.idxs = obj.output.ubloxgps.sc2.trimmed_idxs;
%             obj.gnss.cov = obj.output.ubloxgps.sc2.trimmed_covs;
%             obj.gnss.meas = obj.output.ubloxgps.sc2.trimmed_meas;
        end

        %% Optimization via Solving Non-linear Least Squares
        function obj = optimize(obj)
            disp('-Batch Optimization-')            

            obj.opt.x = lsqnonlin(@obj.cost_func, obj.x0, obj.opt.lb, obj.opt.ub, obj.opt.options);

            [m,n] = size(obj.ins.state_saver);             
            obj.opt.wsf = obj.opt.x(1);            
            obj.opt.states = reshape(obj.opt.x(2:m*n+1),m,n);  
            obj.opt.params = obj.opt.x(m*n+2:end)';

            obj.computeError();

%             obj.opt.info_mat = obj.opt.jac' * obj.opt.jac;
%             
%             % Covariance Extraction using SuiteSparse Toolbox
%             % Need to include for citation in paper
%             
%             disp('-Recovering Covariance using SuiteSparse Toolbox-')
%             [obj.opt.cov, ~] = sparseinv(obj.opt.info_mat); 

            if obj.plot_flag
                obj.plotRes();
            end

%             if strcmp(obj.mode,'full')
%                 % Lane Point reordering
%                 obj.reorderLP();
%         
%             end
        end

        %% Optimization Cost Function
        function res = cost_func(obj, x0)
            wsf = x0(1);
            [m,n] = size(obj.ins.state_saver);                                    
            states = reshape(x0(2:m*n+1),m,n);   
            paramsL = x0(m*n+1+1:m*n+1+4)';
            paramsR = x0(m*n+1+5:end)';

            obj.CreatePrBlock(wsf,states,paramsL,paramsR);
            obj.CreateGNSSBlock(states); 
            obj.CreateMMBlock(states); 
            obj.CreateWSBlock(states,wsf); 
            obj.CreateMEBlock(states,params);
            
            obj.Merge();
            
            res = obj.opt.res;
        end
        
        %% Prior Block and Residual
        function obj = CreatePrBlock(obj,wsf,states,paramsL,paramsR)
            % Initial Vehicle State
            prior_cov = diag([0.1^2 ... % WSS
                              0.01^2 0.01^2 0.01^2 ...
                              0.01^2 0.01^2 0.01^2 ...
                              0.01^2 0.01^2 0.01^2 ...
                              0.3^2 0.3^2 0.3^2 ...
                              0.1^2 0.1^2 0.1^2]);
             
            [m,n] = size(states);
            
%             blk_width = 1 + m*n + length(params);
            
            % Vehicle State Prior
            Pr_vec_veh = zeros(m+1,1);
%             Blk_s = sparse([],[],[],m+1,blk_width);
%             Blk_s(:,1:m+1) = InvMahalanobis(eye(m+1),prior_cov);
            
            Pr_diff = vertcat(-1 + wsf, ...
                              -obj.ins.state_saver(:,1) + states(:,1));
            Pr_vec_veh(1:m+1) = InvMahalanobis(Pr_diff,prior_cov);
            
            % Arc Parameters Prior
            % Only take the initial segment parameters for prior
            init_paramsL = paramsL(1:4);
            prior_cov = diag([0.01^2 0.01^2 0.5^2 1e-4]);
            Pr_vec_arcL = InvMahalanobis(init_paramsL' - obj.ins.paramsL,prior_cov);

            init_paramsR = paramsR(1:4);
            prior_cov = diag([0.01^2 0.01^2 0.5^2 1e-4]);
            Pr_vec_arcR = InvMahalanobis(init_paramsR' - obj.ins.paramsR,prior_cov);

            obj.opt.Pr_vec = [Pr_vec_veh; Pr_vec_arcL; Pr_vec_arcR];
        end
        
        %% GNSS Block and Residual
        function obj = CreateGNSSBlock(obj,states)
            Jx = horzcat(eye(3),zeros(3,12));
            [m,n] = size(states);
            mJ = size(Jx,1);

%             if strcmp(obj.mode,'full')
%                 blk_width = m*n + 2*2*n*obj.prev_num;
%             elseif strcmp(obj.mode,'partial')
%                 blk_width = m*n;
%             end

            blk_height = size(obj.gnss.meas,2)*mJ;
            GNSS_vec = zeros(blk_height,1);
            
%             tmp = m*3;
            d_gnss = size(obj.gnss.meas,2);
            
%             I = zeros(1,tmp*d_gnss); J = I; V = I;

            for i=1:d_gnss
                state_idx = obj.gnss.idxs(i);
                rel_idx = state_idx - obj.dstart+1;
                
                gnss_cov = diag([obj.gnss.cov(1,i), obj.gnss.cov(5,i), obj.gnss.cov(9,i)]);
%                              
%                 [I_s, J_s, V_s] = sparseFormat(mJ*i-mJ+1:mJ*i,m*(rel_idx-1)+1:m*rel_idx,InvMahalanobis(Jx,gnss_cov));
%                 I(tmp*i-(tmp-1):tmp*i) = I_s; J(tmp*i-(tmp-1):tmp*i) = J_s; V(tmp*i-(tmp-1):tmp*i) = V_s;
                
                % Residual
                Z_gnss = obj.gnss.meas(:,i);
                Z_pred = Jx * states(:,rel_idx);
                Z_diff = -Z_gnss + Z_pred;
                
                GNSS_vec(mJ*i-mJ+1:mJ*i) = InvMahalanobis(Z_diff, gnss_cov);
            end
            
%             GNSS_block = sparse(I, J, V, blk_height, blk_width);
%             dummy = sparse([],[],[], blk_height, 1);
            
            obj.opt.GNSS_vec = GNSS_vec;
%             obj.opt.GNSS_block = horzcat(dummy, GNSS_block);

        end

        %% Motion Model Block and Residual
        function obj = CreateMMBlock(obj,states)
            [m,n] = size(states);
            blk_height = m*(n-1);
            MM_vec = zeros(blk_height,1);
            
            ab_std = 0.9e-5;
            wb_std = 2*1e-06;                        

%             if strcmp(obj.mode,'full')
%                 blk_width = m*n + 2*2*n*obj.prev_num;
%             elseif strcmp(obj.mode,'partial')
%                 blk_width = m*n;
%             end
            
%             tmp = 2*m^2;
%             I = zeros(1,tmp*(n-1)); J = I; V = I;

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
                
%                 Jx = JacobianMM(X_curr, u, dt);
%                    
%                 [I_s1, J_s1, V_s1] = sparseFormat(m*i-(m-1):m*i,m*i-(m-1):m*i,InvMahalanobis(Jx, MM_cov));
%                 [I_s2, J_s2, V_s2] = sparseFormat(m*i-(m-1):m*i,m*i+1:m*i+m,InvMahalanobis(-eye(m), MM_cov));
%                 I(tmp*i-(tmp-1):tmp*i) = horzcat(I_s1, I_s2); J(tmp*i-(tmp-1):tmp*i) = horzcat(J_s1, J_s2); V(tmp*i-(tmp-1):tmp*i) = horzcat(V_s1, V_s2);
                        
                X_diff = -X_next + X_pred;
                MM_vec(m*i-(m-1):m*i) = InvMahalanobis(X_diff, MM_cov);
            end
%             MM_block = sparse(I, J, V, blk_height, blk_width);
%             dummy = sparse([],[],[],blk_height,1);
            
            obj.opt.MM_vec = MM_vec;
%             obj.opt.MM_block = horzcat(dummy, MM_block);

        end
        
        %% Wheel Speed Measurement Block and Residual
        function obj = CreateWSBlock(obj,states,wsf)
            m = size(states,1);
            n = size(states,2);
            
            blk_height = 3*n;

%             if strcmp(obj.mode,'full')
%                 blk_width = 1 + m*n + 2*2*n*obj.prev_num;
%             elseif strcmp(obj.mode,'partial')
%                 blk_width = 1 + m*n;
%             end
            
            WS_vec = zeros(blk_height,1);

            cov = 5e-4*eye(3);
%             tmp = 3+3*m;
%             I = zeros(1,tmp*n); J = I; V = I;
            
            for i=1:n
                X = states(:,i);
%                 [Jwsf, Jx] = JacobianWS(X,wsf,obj.lr,obj.lane.xs.w(3,obj.dstart+i-1));
%                 [Iw, Jw, Vw] = sparseFormat(3*i-2:3*i,1,InvMahalanobis(Jwsf,cov));
%                 [Iv, Jv, Vv] = sparseFormat(3*i-2:3*i,1+m*i-(m-1):1+m*i,InvMahalanobis(Jx,cov));
% 
%                 I(tmp*i-(tmp-1):tmp*i) = horzcat(Iw,Iv);
%                 J(tmp*i-(tmp-1):tmp*i) = horzcat(Jw,Jv);
%                 V(tmp*i-(tmp-1):tmp*i) = horzcat(Vw,Vv);

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


%             WS_block = sparse(I,J,V,blk_height,blk_width);

            obj.opt.WS_vec = WS_vec;
%             obj.opt.WS_block = WS_block;
        end

        %% Measurement Model Block and Residual 
        function obj = CreateMEBlock(obj,states,paramsL,paramsR)
            % Only First Segment is implemented, further step is to be done
            % in the future (Multi arc spline measurement model)

            [m,n] = size(states);
            for i=1:n
                X = states(:,i);
                for j=1:obj.prev_num
                    Zpred_l = MeasPred(X,paramsL(1:4),j);
                    Zpred_r = Measpred(X,paramsR(1:4),j);
                end
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
                                  obj.opt.ME_vec, ...
                                  obj.opt.LP_vec);

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
            elseif strcmp(obj.mode,'partial')
                [opt_lml_lon, opt_lml_lat] = lin_to_geo(obj.opt.lml(1,:), obj.opt.lml(2,:), ref);
                [opt_lmr_lon, opt_lmr_lat] = lin_to_geo(obj.opt.lmr(1,:), obj.opt.lmr(2,:), ref);

                p_left = plot(opt_lml_lon, opt_lml_lat,'m--');
                p_right = plot(opt_lmr_lon, opt_lmr_lat,'c--');
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
                legend([p_opt1, p_opt2, p_gt, p_gnss p_left p_right],...
                       'Optimized Trajectory(w GNSS)','Optimized Trajectory(w/o GNSS)',...
                       'Ground Truth','GNSS Measurement(Ublox)',...
                       'Added Left Lane', 'Added Right Lane')
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

            %% Covariance Analysis

%             obj.opt.info_mat = obj.opt.jac' * obj.opt.jac;
%             veh_info = obj.opt.info_mat(1:1+m*n,1:1+m*n);
% 
%             figure(9); hold on; grid on;
%             for i=1:3
%                 idxs = 1+i:m:1+m*n;
%                 covs = sqrt(diag(inv(veh_info(idxs,idxs))));
%                 switch i
%                     case 1
%                         color = 'r';
%                     case 2
%                         color = 'g';
%                     case 3
%                         color = 'b';
%                 end
%                 plot(t, covs, color);
%             end
%             xlabel('Time(s)'); ylabel('m'); title('3D Coordinates Covariance');
%             
%             figure(10);
%             idxs = 10:m:1+m*n;
%             covs = 180/pi*sqrt(diag(inv(veh_info(idxs,idxs))));
%             plot(t, covs);
%             xlabel('Time(s)'); ylabel('Angle(Deg)'); title('Heading Angle Covariance');
%             
%             if strcmp(obj.mode,'full')
%                 lane_info = obj.opt.info_mat(2+m*n:end,2+m*n:end);
%                 
%                 figure(11); hold on; grid on;
%                 for i=1:obj.prev_num
%                     idxs_x = 2*i-1:2*obj.prev_num:2*n*obj.prev_num;
%                     idxs_y = 2*i:2*obj.prev_num:2*n*obj.prev_num;
%                     covs_x = sqrt(diag(inv(lane_info(idxs_x,idxs_x))));
%                     covs_y = sqrt(diag(inv(lane_info(idxs_y,idxs_y))));
%                     
%                     subplot(2,1,1)
%                     plot(t, covs_x)
%                     xlabel('Time(s)'); ylabel('m'); title('Left Lane Point Covariance X');
% 
%                     subplot(2,1,2)
%                     plot(t, covs_y)
%                     xlabel('Time(s)'); ylabel('m'); title('Left Lane Point Covariance Y');
%                 end
%                 
% 
%                 figure(12); hold on; grid on;
%                 for i=1:obj.prev_num
%                     idxs_x = 2*n*obj.prev_num+2*i-1:2*obj.prev_num:4*n*obj.prev_num;
%                     idxs_y = 2*n*obj.prev_num+2*i:2*obj.prev_num:4*n*obj.prev_num;
%                     covs_x = sqrt(diag(inv(lane_info(idxs_x,idxs_x))));
%                     covs_y = sqrt(diag(inv(lane_info(idxs_y,idxs_y))));
%                     
%                     subplot(2,1,1)
%                     plot(t, covs_x)
%                     xlabel('Time(s)'); ylabel('m'); title('Right Lane Point Covariance X');
% 
%                     subplot(2,1,2)
%                     plot(t, covs_y)
%                     xlabel('Time(s)'); ylabel('m'); title('Right Lane Point Covariance Y');
%                 end
%     
%                 
%             end
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
            
%             [obj.opt.c_info, ...
%              obj.opt.reordered_lml_pc, ...
%              obj.opt.reordered_lmr_pc, ...
%              obj.opt.reordered_lml_cov, ...
%              obj.opt.reordered_lmr_cov] = clusterLP(lml_pc, lmr_pc, lml_cov, lmr_cov);
            

            [obj.opt.reordered_lml_pc, obj.opt.reordered_lmr_pc, ...
             obj.opt.reordered_lml_cov, obj.opt.reordered_lmr_cov] = sortLP(lml_pc, lmr_pc, lml_cov, lmr_cov);

            obj.getWeights();
            

%             % Cluster close lane points together
%             obj.opt.cluster_info_l = clusterLP(lml_pc);
%             obj.opt.cluster_info_r = clusterLP(lmr_pc);
%             obj.cluster();
% 
%             obj.opt.reordered_lml_pc = zeros(2,size(obj.opt.cluster_info,2));
%             obj.opt.reordered_lmr_pc = zeros(2,size(obj.opt.cluster_info,2));
% 
%             obj.opt.reordered_lml_cov = zeros(4,size(obj.opt.cluster_info,2));
%             obj.opt.reordered_lmr_cov = zeros(4,size(obj.opt.cluster_info,2));
%                     
%             num = length(obj.opt.cluster_info);
%             for i=1:num
%                 idxs = obj.opt.cluster_info{i};
%                 if length(idxs) == 1 % no close points
%                     obj.opt.reordered_lml_pc(:,i) = lml_pc(:,idxs);
%                     obj.opt.reordered_lml_cov(:,i) = lml_cov(:,idxs);
%                     obj.opt.reordered_lmr_pc(:,i) = lmr_pc(:,idxs);
%                     obj.opt.reordered_lmr_cov(:,i) = lmr_cov(:,idxs);
% 
%                 else % close points exist
%                     lps_l = []; lps_r = [];
%                     covs_l = []; covs_r = [];
% 
%                     m = length(idxs);
% 
%                     for j=1:length(idxs)
%                         lps_l = [lps_l lml_pc(:,idxs(j))];
%                         covs_l = [covs_l lml_cov(:,idxs(j))];
%                         lps_r = [lps_r lmr_pc(:,idxs(j))];
%                         covs_r = [covs_r lmr_cov(:,idxs(j))];
%                     end
%                     
%                     obj.opt.reordered_lml_pc(:,i) = mean(lps_l,2);
%                     obj.opt.reordered_lml_cov(:,i) = mean(covs_l,2)/m;
%                     obj.opt.reordered_lmr_pc(:,i) = mean(lps_r,2);
%                     obj.opt.reordered_lmr_cov(:,i) = mean(covs_r,2)/m;
%                 end
%             end
        end

        %% Reshape x0 format for Lane Point
        function res = InvreshapeLP(obj,InvLP)
            n = size(obj.ins.state_saver,2);            
            res = zeros(size(obj.ins.lml_saver));
            
            tmp = 2*obj.prev_num;
            
            for j=1:n
                res(2*j-1:2*j,:) = reshape(InvLP(tmp*j-(tmp-1):tmp*j),2,obj.prev_num);
            end            
        end
        
        %% Compute Vehicle Frame Covariance
        function Cov = getCov(obj, lat_std)
            Cov = zeros(4*(obj.dend-obj.dstart+1), obj.prev_num); % Exclude the initial point
            for i=obj.dstart:obj.dend
                rel_idx = i - obj.dstart+1;
                for j=1:obj.prev_num
                    Cov(4*rel_idx-3:4*rel_idx,j) = reshape([0.03^2 0; 0 lat_std(i,j)^2],4,1);
                end
            end
        end
        
        %% Cluster Indices of Left & Right Lanes 
        function obj = cluster(obj)
            % Reference: Left Lane
            search_idx_l = 1;
            search_idx_r = 1;
            obj.opt.cluster_info = {};
            obj.opt.left_accum = [];
            obj.opt.right_accum = [];

            while search_idx_l <= length(obj.opt.cluster_info_l)
                left_idxs = obj.opt.cluster_info_l{search_idx_l};
                right_idxs = obj.opt.cluster_info_r{search_idx_r};
                
                obj.opt.left_accum = [obj.opt.left_accum left_idxs];
                obj.opt.right_accum = [obj.opt.right_accum right_idxs];
                
                L_l = length(obj.opt.left_accum); L_r = length(obj.opt.right_accum);
                
                if length(left_idxs) > length(right_idxs)
                    rep_idxs = right_idxs;
                else
                    rep_idxs = left_idxs;
                end

                while L_l~=L_r
                    if L_l > L_r 
                        search_idx_r = search_idx_r + 1;
                        right_idxs = obj.opt.cluster_info_r{search_idx_r};
                        obj.opt.right_accum = [obj.opt.right_accum right_idxs];
                        L_r = length(obj.opt.right_accum);
                        rep_idxs = [rep_idxs right_idxs];
                    else
                        search_idx_l = search_idx_l + 1;
                        left_idxs = obj.opt.cluster_info_l{search_idx_l};
                        obj.opt.left_accum = [obj.opt.left_accum left_idxs];
                        L_l = length(obj.opt.left_accum);
                        rep_idxs = [rep_idxs left_idxs];
                    end
                end
                
                obj.opt.cluster_info = [obj.opt.cluster_info {rep_idxs}];

                search_idx_l = search_idx_l + 1;
                search_idx_r = search_idx_r + 1;
            end
        end
        
        %% Convert Covariance Matrix to weights
        function obj = getWeights(obj)
            n = size(obj.opt.reordered_lml_cov,2);
            obj.opt.w_l = zeros(1,n);
            obj.opt.w_r = zeros(1,n);
            p = 0.95;
            
            for i=1:n

                cov_l = reshape(obj.opt.reordered_lml_cov(:,i),2,2);
                cov_r = reshape(obj.opt.reordered_lmr_cov(:,i),2,2);
                eig_l = min(eig(cov_l));
                eig_r = min(eig(cov_r));
                
                obj.opt.w_l(i) = chi2inv(p,2) * sqrt(eig_l);
                obj.opt.w_r(i) = chi2inv(p,2) * sqrt(eig_r);
            end
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

%% Measurement Prediction
function Zpred = MeasPred(X,params,idx)
    x0 = params(1); y0 = params(2); tau0 = params(3); kappa = params(4);
    xc = x0 - 1/kappa * sin(tau0);
    yc = y0 + 1/kappa * cos(tau0);

    xv_org = X(1); yv_org = X(2); psi = X(9);
    xv = xv_org + 10*(idx-1)*cos(psi);
    yv = yv_org + 10*(idx-1)*sin(psi);

    delx = xc - xv; dely = yc - yv;
    xcb = cos(psi) * delx + sin(psi) * dely;
    ycb = -sin(psi) * delx + cos(psi) * dely;

    if ycb < 0
        Zpred = ycb + sqrt(1/kappa^2 - xcb^2);
    else
        Zpred = ycb - sqrt(1/kappa^2 - xcb^2);
    end
end

%% Lane Point Prediction
function res = LanePred(Xk, Lij, Lijp1, Lipkj, idx)
    xk = Xk(1); yk = Xk(2); psik = Xk(9);
    lijx = Lij(1); lijy = Lij(2);
    lijp1x = Lijp1(1); lijp1y = Lijp1(2);
    lipkjx = Lipkj(1); lipkjy = Lipkj(2);

    delx_ij = lijx - xk; delx_ijp1 = lijp1x - xk; delx_ipkj = lipkjx - xk;
    dely_ij = lijy - yk; dely_ijp1 = lijp1y - yk; dely_ipkj = lipkjy - yk;
    rot = [cos(psik) sin(psik); -sin(psik) cos(psik)];
    
    Lijb = rot * [delx_ij; dely_ij];
    Lijp1b = rot * [delx_ijp1; dely_ijp1];
    Lipkjb = rot * [delx_ipkj; dely_ipkj];
    
    x_pred = 10*(idx-1);
    y_pred = (Lijp1b(2) - Lijb(2))/(Lijp1b(1) - Lijb(1)) * (x_pred - Lijb(1)) + Lijb(2);
    L_pred = [x_pred; y_pred];
    
    res = L_pred - Lipkjb;
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

%% Rotate Covariance Matrix
function CovMat = rotCov(Cov, alpha)
    
    n = size(Cov,2);
    CovMat = zeros(size(Cov));
    Rot = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];

    for i=1:n
        cov = reshape(Cov(:,i),2,2);
        covrot = Rot * cov * Rot';
        CovMat(:,i) = reshape(covrot,4,1);
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

%% Cluster close points
function [c_info, re_lml, re_lmr, re_cov_l, re_cov_r] = clusterLP(LP_l, LP_r, cov_l, cov_r)
    LP_cpyl = LP_l; LP_cpyr = LP_r;
    LP_cpyl(3,:) = 1:size(LP_l,2); % add absolute index at the last row
    LP_cpyr(3,:) = 1:size(LP_r,2);
    
    lp_l = LP_cpyl(1:2,1); % initialized point for clustering
    lp_r = LP_cpyr(1:2,1);
    abs_idx = 1; % absolute index
    LP_cpyl = LP_cpyl(:,2:end);
    LP_cpyr = LP_cpyr(:,2:end);
    thres = 0.8;
    
    re_cov_l = [];
    re_cov_r = [];
    re_lml = [];
    re_lmr = [];
    c_info = {};

    while size(LP_cpyl,2) > 0
        delx_l = LP_cpyl(1,:) - lp_l(1);
        dely_l = LP_cpyl(2,:) - lp_l(2);
        delx_r = LP_cpyr(1,:) - lp_r(1);
        dely_r = LP_cpyr(2,:) - lp_r(2);
        
        d = sqrt(delx_l.^2 + dely_l.^2) + sqrt(delx_r.^2 + dely_r.^2);        
        [d_min,new_idx] = min(d);
             
        if d_min > thres % No close lane point
            c_info = [c_info, {abs_idx}];
            re_cov_l = [re_cov_l cov_l(:,abs_idx)];
            re_cov_r = [re_cov_r cov_r(:,abs_idx)];
            re_lml = [re_lml LP_l(:,abs_idx)];
            re_lmr = [re_lmr LP_r(:,abs_idx)];
            
            lp_l = LP_cpyl(1:2,new_idx);
            lp_r = LP_cpyr(1:2,new_idx);
            abs_idx = LP_cpyl(3,new_idx);
            LP_cpyl(:,new_idx) = [];
            LP_cpyr(:,new_idx) = [];
            
        else % Close lane point exists to current point
            idxs = find(d<thres);
            
            abs_idxs = abs_idx;
            lps_l = lp_l;
            lps_r = lp_r;
            covs_l = cov_l(:,abs_idx);
            covs_r = cov_r(:,abs_idx);
            for i=1:length(idxs)
                abs_idxs = [abs_idxs LP_cpyl(3,idxs(i))];
                lps_l = [lps_l LP_cpyl(1:2,idxs(i))];
                lps_r = [lps_r LP_cpyr(1:2,idxs(i))];
                covs_l = [covs_l cov_l(:,LP_cpyl(3,idxs(i)))];
                covs_r = [covs_r cov_r(:,LP_cpyl(3,idxs(i)))];
            end

            c_info = [c_info, {abs_idxs}];
            
            lp_l = mean(lps_l,2);
            lp_r = mean(lps_r,2);
            covs_l = 1/length(idxs) * mean(covs_l,2);
            covs_r = 1/length(idxs) * mean(covs_r,2);
            
            re_cov_l = [re_cov_l covs_l];
            re_cov_r = [re_cov_r covs_r];
            re_lml = [re_lml lp_l];
            re_lmr = [re_lmr lp_r];
            
            LP_cpyl = LP_cpyl(:,d>=thres);
            LP_cpyr = LP_cpyr(:,d>=thres);
            
            delx_l = LP_cpyl(1,:) - lp_l(1);
            dely_l = LP_cpyl(2,:) - lp_l(2);
            delx_r = LP_cpyr(1,:) - lp_r(1);
            dely_r = LP_cpyr(2,:) - lp_r(2);
            d = sqrt(delx_l.^2 + dely_l.^2) + sqrt(delx_r.^2 + dely_r.^2);        
            [~,new_idx] = min(d);
            
            lp_l = LP_cpyl(1:2,new_idx);
            lp_r = LP_cpyr(1:2,new_idx);
            abs_idx = LP_cpyl(3,new_idx);
            LP_cpyl(:,new_idx) = [];
            LP_cpyr(:,new_idx) = [];
        end  
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
end

%% Get circle from 3 points
function kappa = getCircle(LPs)
    if size(LPs,2) ~= 3
        error('Input should be 2 * 3 Array Form (3 Points)')
    end
    
    lp1 = LPs(:,1); x1 = lp1(1); y1 = lp1(2);
    lp2 = LPs(:,2); x2 = lp2(1); y2 = lp2(2);
    lp3 = LPs(:,3); x3 = lp3(1); y3 = lp3(2);
    
    m12 = (y2 - y1)/(x2 - x1);
    m23 = (y3 - y2)/(x3 - x2);
    
    lp12 = 1/2 * [x1+x2;y1+y2];
    lp23 = 1/2 * [x2+x3;y2+y3];

    xc = (lp12(1)/m12 + lp12(2) - lp23(1)/m23 - lp23(2))/(1/m12 - 1/m23);
    yc = -1/m12 * (xc - lp12(1)) + lp12(2);
    diff = [xc; yc] - lp1;
    R = sqrt(diff' * diff);
    vec12 = lp2 - lp1; vec23 = lp3 - lp2;
    th12 = atan2(vec12(2),vec12(1)); th23 = atan2(vec23(2),vec23(1));
    if th23 - th12 >= 0
        kappa = 1/R;
    else
        kappa = -1/R;
    end
end
