classdef FootStepGenerator < handle

    methods (Access = public)

        function obj = FootStepGenerator(input, sim_parameters)
            
            % constructor
            obj.stride_length_x = 0.2;
            obj.lateral_displacement_y = 0.2;
            obj.input = input;
            obj.P_ = obj.Tp_ / obj.delta_;
            obj.plan_ = obj.input.footstep_plan;
            obj.SupportFoot = obj.plan_.starting_sf;
            obj.v_inputs_ = obj.input.footstep_plan.input_velocity;
            ref_pose = obj.input.footstep_plan.positions(2:end,:);
            obj.P_traj =[];
            
            for k = 1: length(obj.v_inputs_)
                P_traj_ = IntegrateProfile(obj, k);
                obj.P_traj = [obj.P_traj; P_traj_];
            end
            
            obj.t_steps_inputs_ = obj.input.footstep_plan.timings(2:end);
            if(~isempty(obj.t_steps_inputs_))
                obj.t_steps_inputs_(1) = max(obj.Ts_min_, obj.t_steps_inputs_(1));
            end
            
            for k = 1 : length(ref_pose)
                if(abs(ref_pose(k,3)) > pi)
                    ref_pose(k,3) = ref_pose(k,3) - (ref_pose(k,3) / abs(ref_pose(k,3))) * 2 * pi ;
                end
            end
            obj.pose_reference_ = ref_pose;
            
            obj.sim_parameters = sim_parameters;
            obj.options = optimoptions(@quadprog, 'Display','final-detailed');
        end
        
        function GetStepsTiming(obj)
            obj.StepTimings_ = obj.t_steps_inputs_;
            obj.StepTimings_idx = [];
            obj.FootSteps_idx = [];
            if(~isempty(obj.pose_reference_))
                P_f_im1 = IntegrateProfile(obj,1);
                P_f_i = obj.pose_reference_(1,:);
                d_step = norm([obj.d_h_x, obj.d_h_y])/2;
                nsteps = floor(norm(P_f_im1(1:2) - P_f_i(1:2))/d_step + 1);
                if(nsteps < length(obj.StepTimings_))
                    obj.P_traj = GetRefTrajectory(obj, P_f_im1, P_f_i, obj.StepTimings_(nsteps));
                    while (length(obj.P_traj) < obj.P_)
                        obj.P_traj = [obj.P_traj; obj.P_traj(end,:)];
                    end
                else
                    obj.P_traj = GetRefTrajectory(obj, P_f_im1, P_f_i, obj.Tp_);
                end
            end
            for i = 1: length(obj.StepTimings_)
                t_i = obj.StepTimings_(i);
                obj.StepTimings_idx = [obj.StepTimings_idx floor(t_i/obj.delta_)];
                obj.FootSteps_idx = [obj.FootSteps_idx -1];
            end
            while(obj.StepTimings_(length(obj.StepTimings_) - 1) > obj.Tp_ - obj.delta_)
                obj.StepTimings_ = obj.StepTimings_(1 : length(obj.StepTimings_) - 1);
                obj.FootSteps_idx = obj.FootSteps_idx(1 : length(obj.FootSteps_idx) - 1);
                obj.StepTimings_idx = obj.StepTimings_idx(1 : length(obj.StepTimings_idx) - 1);
            end
            obj.F_ = length(obj.StepTimings_);

        end
        
        function  compute_plan(obj)
            obj.plan_.positions = [];
            GetStepsTiming(obj);
            obj.Theta_f_ = zeros(obj.F_, 1);
            if (obj.F_ <= 0)
                error("[footsteps_planner::compute_plan] No step to compute");
                obj.plan_.positions = [];
                return ;
            end
            
            %% Solving QP 1 for orientation 
            
            Dtheta_upper_lim = ones(obj.F_,1) * obj.max_theta;
            Dtheta_lower_lim = ones(obj.F_,1) * -obj.max_theta;
            Delta = eye(obj.F_);
            obj.A_eq = zeros(obj.F_);
            obj.b_eq = zeros(obj.F_,1);
            for k = 1: obj.F_
                if (k == 1)
                    Dtheta_upper_lim(k) = Dtheta_upper_lim(k) + obj.plan_.sf_pos(3) ;
                    Dtheta_lower_lim(k) = Dtheta_lower_lim(k) + obj.plan_.sf_pos(3) ;
                else
                    Delta(k, k -1) = -1;
                end
                if(obj.FootSteps_idx(k) ~= -1)
                      theta_cstr = obj.pose_reference_(obj.FootSteps_idx(k), 3);
                      
                      if(obj.P_traj(obj.StepTimings_idx(k),3) - obj.pose_reference_(obj.FootSteps_idx(k),3) > pi)
                          
                        theta_cstr = theta_cstr - obj.pose_reference_(obj.FootSteps_idx(k),3) / abs(obj.pose_reference_(obj.FootSteps_idx(k),3)) * 2 * pi;
                      end
                      obj.A_eq(k, k) = 1;
                      obj.b_eq(k) = theta_cstr;
                end
            end
            obj.A_ineq = zeros(2 * obj.F_, obj.F_);
            obj.b_ineq = zeros(2 * obj.F_, 1);
            
            obj.A_ineq = [Delta ; -Delta];
            obj.b_ineq = [Dtheta_upper_lim; -Dtheta_lower_lim];
            
            b = zeros(obj.F_,1);
            for i = 1 : obj.F_
                b(i) = obj.P_traj(obj.StepTimings_idx(i),3);
            end
            obj.H = eye(obj.F_);
            obj.f = -b;
            theta = solveQP(obj);
            if (isempty(theta))
                error("[footsteps planner] Step theta failed");
                 obj.plan_.positions = [];
                return ;
            end
            for k = 1 : length(theta)
                if abs(theta(k)) > pi
                    theta_f_(k,1) = theta(k) - (theta(k) / abs(theta(k))) * 2 * pi;
                else
                    theta_f_(k,1) = theta(k);
                end
            end            
            %% Solving QP 2 for placement
            Delta = eye(2 * obj.F_);
            obj.A_eq = zeros(2 * obj.F_);
            obj.b_eq = zeros(2 * obj.F_, 1);
            
            sgn = -1.0;
            if (obj.SupportFoot == "right")
                sgn = 1.0;
            end
            l = obj.l_;
            if (obj.d_h_y /2 > obj.l_ - 0.5 * obj.d_min)
                l = 0.5 * (obj.d_h_y + obj.d_min);
            end
            
            Normal_Vec = [];
            cstr_vec = [];
            
            theta_0 = obj.P_traj(1,3);
            R_theta_0 = rotz(-rad2deg(theta_0));
            R_theta_0 = R_theta_0(1:2,1:2);
            
            kinematic_rectangle = Admissible_Region([0,0,obj.plan_.sf_pos(3)],[obj.d_h_x, obj.d_h_y, 0]);
            bcstr = kinematic_rectangle.offset + kinematic_rectangle.polygone_Normals'...
                  * (obj.plan_.sf_pos(1:2) + sgn * R_theta_0 * [0;l]);
            
            Normal_Vec = [Normal_Vec; kinematic_rectangle.polygone_Normals'];
            cstr_vec = [cstr_vec; bcstr];
            
            for k = 2 : obj.F_
                theta_k = obj.P_traj(obj.StepTimings_idx(k),3);
                theta_f_km1 = theta_f_(k-1);

                kinematic_rectangle = Admissible_Region([0, 0, theta_f_km1], [obj.d_h_x, obj.d_h_y, 0]);
                
                R = rotz(rad2deg(-theta_f_km1));
                R_Theta_f_km1 = R(1:2, 1:2);
                bcstr = kinematic_rectangle.offset + sgn * (1- 2 * mod(k-1, 2)) * kinematic_rectangle.polygone_Normals' ...
                        * R_Theta_f_km1 * [0; l];

                Normal_Vec = [Normal_Vec; kinematic_rectangle.polygone_Normals'];
                cstr_vec = [cstr_vec; bcstr];
            end

            if (obj.FootSteps_idx(1) ~= -1) %%%%%%% revise %%%%%%ù
              obj.A_eq(1,1) = 1;
              obj.A_eq(2,2) = 1;
              obj.b_eq(1) = obj.plan_.positions(obj.FootSteps_idx(1),1);
              obj.b_eq(2) = obj.plan_.positions(obj.FootSteps_idx(1),2);
            end
            
            for k = 2 : obj.F_
                Delta(2 * k - 1, 2 * (k-1) - 1) = -1;
                Delta(2 * k    , 2 *  k - 2) = -1;
                
                if (obj.FootSteps_idx(k) ~= -1) %%%%%%% revise %%%%%%ù
                  obj.A_eq(2 * k + 1, 2 * k + 1) = 1;
                  obj.A_eq(2 * k + 2, 2 * k + 2) = 1;
                  obj.b_eq(2 * k + 1) = obj.plan_.positions(obj.FootSteps_idx(k + 1),1);
                  obj.b_eq(2 * k + 2) = obj.plan_.positions(obj.FootSteps_idx(k + 1),2);
                end
            end
           
            
            N_cstr = length(Normal_Vec);
            
            obj.A_ineq = zeros(N_cstr, 2 *obj.F_);
            obj.b_ineq = zeros(N_cstr,1);
            
            step_indx = 1;
            cstr_indx = 1;
            
            for i_ineq = 1: length(Normal_Vec) / 4 
                n_vec = Normal_Vec(4 * (i_ineq - 1) + 1: 4 * i_ineq, :);
                ineq = cstr_vec(4 * (i_ineq - 1) + 1: 4 * i_ineq);
            
                for cstr = 1: length(n_vec)
                    obj.A_ineq(cstr_indx + cstr, step_indx) = n_vec(cstr, 1);
                    obj.A_ineq(cstr_indx + cstr, step_indx + 1) = n_vec(cstr, 2);

                    obj.b_ineq(cstr_indx + cstr) = ineq(cstr);
                end
                step_indx = step_indx + 2;
                cstr_indx = cstr_indx + length(n_vec);
            end
            obj.A_ineq = obj.A_ineq * Delta;
            b = zeros(2 * obj.F_,1);
            for i = 1: obj.F_-1
                theta_i = obj.P_traj(i,3);
                R = rotz(rad2deg(-theta_i));
                R_theta_i = R(1:2,1:2);
                temp = (1 - 2 * (mod(i-1,2)));
                dl = (sgn * temp * R_theta_i * [0; obj.l_/2]);
                b(2 * i - 1) = obj.P_traj(obj.StepTimings_idx(i),1) + dl(1);
                b(2 * i) = obj.P_traj(obj.StepTimings_idx(i),2) + dl(2);
            end
            obj.H = eye(2 * obj.F_);
            obj.f = -b;
            
            XY = solveQP(obj);
            if(isempty(XY))
                error("[footsteps planner] Step QP failed");
                obj.plan_.positions = [];
                return ;
            else
                for k = 1: obj.F_
                    xf = XY(2 * k  - 1);
                    yf = XY(2 * k);
                    obj.plan_.positions(k,:) = [xf, yf, theta_f_(k)];
                end
            end
        end
        
        function output = GetRefTrajectory(obj, P_s_0, P_s_1, duration)
             output = [];
             R = rotz(rad2deg(P_s_0(3)));
             R_s0_0 = R(1:2,1:2)';
             shuffle = false;
             backward = false;
             temp = abs(R_s0_0' * (P_s_1(1:2)' - P_s_0(1:2)'));
             if (temp(1) < 3e-1)
                 shuffle = true;
             end
             temp = R_s0_0 * (P_s_1(1:2)' - P_s_0(1:2)');
             if (temp(1) < 0  && ~shuffle)
                 backward = true;
             end
             init_pos = P_s_0(1:2)';
             target_pos = P_s_1(1:2)';
             if(~shuffle)
                 init_ori = [cos(P_s_0(3)), sin(P_s_0(3))];
                 target_ori = [cos(P_s_1(3)), sin(P_s_1(3))];
                 if(backward)
                     init_ori = [-cos(P_s_0(3)), -sin(P_s_0(3))];
                     target_ori = [-cos(P_s_1(3)), -sin(P_s_1(3))];
                 end
                 PolynomeCoeff(obj, init_pos, init_ori, target_pos, target_ori);
             end
             N = floor(duration / obj.delta_);
             for k= 0 :N
                 t = k / N;
                 if(~shuffle)
                     path = computePath(obj, t);
                     pos_t = path(1:2);
                     ori_t = path(3:4);
                     theta = atan2(ori_t(2), ori_t(1));
                     if(norm(init_pos - target_pos)> 5e-2)
                        theta = P_s_1(3);
                     end
                     if(backward)
                        sgn = theta / abs(theta);
                        output = [output; [pos_t, mod(theta - sgn * pi, 2 * pi)]];
                     else
                        output = [output; [pos_t, mod(theta, 2 * pi)]];
                     end
                 else
                     output = [output; [P_s_0(1:2) + t * (P_s_1(1:2) - P_s_0(1:2)), P_s_1(3)]];
                 end
             end
        end
        
        function  solution = solveQP(obj)
          
           solution = quadprog(obj.H, obj.f, obj.A_ineq, obj.b_ineq, obj.A_eq, obj.b_eq , [], [], [], obj.options);
            
        end
        
        function output = IntegrateProfile(obj, k_end)
            k_start = 1;
            N = size(obj.P_traj);
            if(k_end - 1 < N(1))
                output = obj.P_traj(k_end,:);
                return;
            end
            output = obj.plan_.positions(1,:)';
            if (isempty(obj.P_traj))
                R = rotz(rad2deg(-output(3)));                
                if(obj.SupportFoot == "right")
                    output(1:2) = output(1:2) + R(1:2,1:2) * [0; obj.l_/ 2];
                else
                    output(1:2) = output(1:2) - R(1:2,1:2) * [0; obj.l_/ 2];
                end
            else
                endPtraj = obj.P_traj(end,:);
                output = endPtraj';
                k_start = N(1);
            end
            for i = k_start:min(k_end, length(obj.v_inputs_))
               R = rotz(rad2deg(-output(3)));
               output(3) = output(3) + obj.v_inputs_(i, 3) * obj.delta_;
               output(1:2) = output(1:2) + R(1:2,1:2) * obj.v_inputs_(i, 1:2)' * obj.delta_;
            end
               output = output';
        end

        function R = RotZ(angle)
           R = [cos(angle), -sin(angle),  0;...
                sin(angle),  cos(angle),  0; ...
                         0,           0,  1];
        end
        
        function PolynomeCoeff(obj,init_pos, initVel, target_pos, targetVel)
            obj.c0 = init_pos;
            obj.c1 = initVel;
            obj.c2 = 3 * (target_pos - init_pos) - 2 * initVel - targetVel;
            obj.c3 = -2 * (target_pos - init_pos) + initVel + targetVel;
        end
        
        function path = computePath(obj,t)
            pos = obj.c0 + t * (obj.c1 + t * (obj.c2 + t * obj.c3));
            vel = obj.c1 + t * (2 * obj.c2 + 3 * t * obj.c3);
            tangent = vel / norm(vel);
            path = [pos, tangent];
        end
        
        function plan = GetFootStepPlan(obj)
            plan = obj.plan_.positions;
        end
        
        function timing = GetFootStepTiming(obj)
            timing = obj.StepTimings_;
        end
        
    end

    properties (Access = private)
       Ts_min_ = 0.8; % Step Time lenght limits
       Ts_max_ = 2; % Step Time lenght limits
       l_ = 0.2; % Distance between foot
       d_min = 0.34; % min distance between the foot center and the kinematic constraint;
       Tp_ = 10; % Preview horizon time
       delta_ = 0.1; % t_k - t_k-1
       d_h_x = 0.2; % Next step tolerance zone
       d_h_y = 0.05; % Next step tolerance zone
       v_ = 0.1; % Cruise Parameters
       max_theta = 3.14 / 6; % Max angle between two steps
       P_ = 100; % Preview horizon time indexes
       F_ = 1; % footsteps number
       Ts_ = 5.0; % Cruise Parameters
       robot_height_ = 150; % in cm
       theta_offset_ = 0;
        
        % building blocks
        input;
        sim_parameters;
        state;
        % variables
        H;
        f;
        delta;
        A_ineq;
        b_ineq;
        A_eq;
        b_eq;
        options;
        P_traj;
        stride_length_x;
        lateral_displacement_y;
        duration;
        P_s_0;
        P_s_1;
        c0; c1; c2; c3;
        v_inputs_;
        pose_reference_;
        t_steps_inputs_;
        Theta_f_;
        plan_;
        StepTimings_idx;
        FootSteps_idx;
        SupportFoot;
        StepTimings_
        % buffers

        
    end
    
end
