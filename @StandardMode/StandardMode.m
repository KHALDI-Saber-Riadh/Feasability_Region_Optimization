classdef StandardMode < FeasibilityDrivenBase & handle
    
    methods (Access = public)
        
        function obj = StandardMode(input)
            
            obj.input = input;
            obj.x_u_m = 0;
            obj.x_u_M = 0;
            obj.y_u_m = 0;
            obj.y_u_M = 0;
            obj.feasibility_region = zeros(8, obj.input.kar.number_of_subregions);
            obj.feasibility_region(:, 1) = [obj.x_u_m; obj.x_u_M; obj.y_u_m; obj.y_u_M; 0; 0; 0; 0];
            obj.centerline_multiplier = obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta * ...
                                        exp( - obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta * (0 : obj.input.scheme_parameters.C - 1));
            obj.tail_multiplier = obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta * ...
                                        exp( - obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta * (obj.input.scheme_parameters.C : obj.input.scheme_parameters.P - 1));            
            
            obj.zmp_tracking_weight = 10;
            obj.P_matrix = obj.input.scheme_parameters.delta * tril( ones(obj.input.scheme_parameters.C, obj.input.scheme_parameters.C) );
            
            obj.H = eye(obj.input.scheme_parameters.C, obj.input.scheme_parameters.C) + ...
                    obj.zmp_tracking_weight * obj.P_matrix' * obj.P_matrix;
            obj.f = zeros(obj.input.scheme_parameters.C, 1);
            obj.A_ineq = zeros(obj.input.scheme_parameters.C * 2, obj.input.scheme_parameters.C);
            obj.A_ineq(1:obj.input.scheme_parameters.C, 1:obj.input.scheme_parameters.C) = obj.input.scheme_parameters.delta * tril( ones(obj.input.scheme_parameters.C, obj.input.scheme_parameters.C));
            obj.A_ineq(obj.input.scheme_parameters.C + 1 : 2 * obj.input.scheme_parameters.C, 1:obj.input.scheme_parameters.C) = - obj.input.scheme_parameters.delta * tril( ones(obj.input.scheme_parameters.C, obj.input.scheme_parameters.C));
            obj.b_ineq = zeros(obj.input.scheme_parameters.C * 2, 1);
            obj.A_eq = zeros(1, obj.input.scheme_parameters.C);
            for i = 1 : obj.input.scheme_parameters.C
                 obj.A_eq(1, i) = (1 / obj.input.scheme_parameters.eta) * ...
                                  ( exp( - (i-1) * obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta) - ...
                                    exp( - i * obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta) );
            end
            obj.b_eq = zeros(1);
            
            obj.restriction_builder = BuildRestrictionFunction(obj.input);
            obj.restriction_x = obj.restriction_builder.getRestrictionX();
            obj.restriction_y = obj.restriction_builder.getRestrictionY();
            
            obj.options = optimoptions(@quadprog, 'Display','off');
            
        end
        
        function [u, ftstp] = solve(obj, state, input)
            
            obj.input = input;    
            
            % x component           
            obj.b_ineq(1:obj.input.scheme_parameters.C, 1) = - state.x(3,1) + obj.input.footstep_plan.zmp_centerline_x + obj.input.scheme_parameters.d_zxf - obj.restriction_x;
            obj.b_ineq(obj.input.scheme_parameters.C : 2 * obj.input.scheme_parameters.C - 1, 1) = + state.x(3,1) - obj.input.footstep_plan.zmp_centerline_x + obj.input.scheme_parameters.d_zxb - obj.restriction_x;
            
            obj.b_eq = state.x(1,1) + state.x(2,1) / obj.input.scheme_parameters.eta ...
                       - state.x(3,1) * (1 - exp( - obj.input.scheme_parameters.eta * obj.input.scheme_parameters.T_c)) ...
                       - obj.tail_multiplier * obj.input.footstep_plan.tail_x ...
                       - obj.input.footstep_plan.tail_x(end,1) * exp( - obj.input.scheme_parameters.eta * obj.input.scheme_parameters.T_p) ...
                       + state.w_bar(1,1) / obj.input.scheme_parameters.eta ^ 2;
        
            obj.f = - 2 * obj.zmp_tracking_weight * (obj.input.footstep_plan.zmp_centerline_x' - state.x(3,1)) * obj.P_matrix;
            solution = quadprog(obj.H, obj.f, obj.A_ineq, obj.b_ineq, obj.A_eq, obj.b_eq, [], [], [], obj.options);
            u(1,1) = solution(1);
            
            % y component
            obj.b_ineq(1:obj.input.scheme_parameters.C, 1) = - state.y(3,1) + obj.input.footstep_plan.zmp_centerline_y + obj.input.scheme_parameters.d_zy / 2 - obj.restriction_y;
            obj.b_ineq(obj.input.scheme_parameters.C : 2 * obj.input.scheme_parameters.C - 1, 1) = + state.y(3,1) - obj.input.footstep_plan.zmp_centerline_y + obj.input.scheme_parameters.d_zy / 2 - obj.restriction_y;
            
            obj.b_eq = state.y(1,1) + state.y(2,1) / obj.input.scheme_parameters.eta ...
                       - state.y(3,1) * (1 - exp( - obj.input.scheme_parameters.eta * obj.input.scheme_parameters.T_c)) ...
                       - obj.tail_multiplier * obj.input.footstep_plan.tail_y ...
                       - obj.input.footstep_plan.tail_y(end,1) * exp( - obj.input.scheme_parameters.eta * obj.input.scheme_parameters.T_p) ...
                       + state.w_bar(2,1) / obj.input.scheme_parameters.eta ^ 2;
          
            obj.f = - 2 * obj.zmp_tracking_weight * (obj.input.footstep_plan.zmp_centerline_y' - state.y(3,1)) * obj.P_matrix;
            solution = quadprog(obj.H, obj.f, obj.A_ineq, obj.b_ineq, obj.A_eq, obj.b_eq, [], [], [], obj.options);
            u(2,1) = solution(1);
            
            % next footstep is given by the plan
            obj.index = state.world_time_iter - round( obj.input.footstep_plan.timings(state.footstep_counter, 1) / obj.input.scheme_parameters.delta ) + 1;
            if obj.index <= obj.input.footstep_plan.ds_samples 
                ftstp = obj.input.footstep_plan.positions(state.footstep_counter + 1, 1:3)';
            else
                ftstp = obj.input.footstep_plan.positions(state.footstep_counter + 2, 1:3)';
            end
 
        end
        
        function is_feasible = feasibilityCheck(obj, state, input)
          
            obj.input = input;
            obj.feasibility_region = zeros(8, obj.input.kar.number_of_subregions);
            obj.x_u_m = obj.centerline_multiplier * (obj.input.footstep_plan.zmp_centerline_x - obj.input.scheme_parameters.d_zxb + obj.restriction_x) ...
                        + obj.tail_multiplier * obj.input.footstep_plan.tail_x ...
                        - state.w_bar(1,1) / obj.input.scheme_parameters.eta ^ 2;
            obj.x_u_M = obj.centerline_multiplier * (obj.input.footstep_plan.zmp_centerline_x + obj.input.scheme_parameters.d_zxf - obj.restriction_x) ...
                        + obj.tail_multiplier * obj.input.footstep_plan.tail_x ...
                        - state.w_bar(1,1) / obj.input.scheme_parameters.eta ^ 2; 
            obj.y_u_m = obj.centerline_multiplier * (obj.input.footstep_plan.zmp_centerline_y - obj.input.scheme_parameters.d_zy/2 + obj.restriction_y) ...
                        + obj.tail_multiplier * obj.input.footstep_plan.tail_y ...
                        - state.w_bar(2,1) / obj.input.scheme_parameters.eta ^ 2; 
            obj.y_u_M = obj.centerline_multiplier * (obj.input.footstep_plan.zmp_centerline_y + obj.input.scheme_parameters.d_zy/2 - obj.restriction_y) ...
                        + obj.tail_multiplier * obj.input.footstep_plan.tail_y ...
                        - state.w_bar(2,1) / obj.input.scheme_parameters.eta ^ 2;  
            obj.feasibility_region(:, 1) = [obj.x_u_m; obj.x_u_M; obj.y_u_m; obj.y_u_M; 0; 0; 0; 0];
            
            is_feasible = false;
            obj.x_u = state.x(1,1) + state.x(2,1) / obj.input.scheme_parameters.eta;
            obj.y_u = state.y(1,1) + state.y(2,1) / obj.input.scheme_parameters.eta; 
            if (obj.x_u >= obj.x_u_m && obj.x_u <= obj.x_u_M) && (obj.y_u >= obj.y_u_m && obj.y_u <= obj.y_u_M)
                is_feasible = true;    
            end
            
        end
        
        function obj = computeFeasibilityRegion(obj, state, input)
            
            obj.input = input;
            obj.feasibility_region = zeros(8, obj.input.kar.number_of_subregions);
            obj.x_u_m = obj.centerline_multiplier * (obj.input.footstep_plan.zmp_centerline_x - obj.input.scheme_parameters.d_zxb + obj.restriction_x) ...
                        + obj.tail_multiplier * obj.input.footstep_plan.tail_x ...
                        - state.w_bar(1,1) / obj.input.scheme_parameters.eta ^ 2;
            obj.x_u_M = obj.centerline_multiplier * (obj.input.footstep_plan.zmp_centerline_x + obj.input.scheme_parameters.d_zxf - obj.restriction_x) ...
                        + obj.tail_multiplier * obj.input.footstep_plan.tail_x ...
                        - state.w_bar(1,1) / obj.input.scheme_parameters.eta ^ 2; 
            obj.y_u_m = obj.centerline_multiplier * (obj.input.footstep_plan.zmp_centerline_y - obj.input.scheme_parameters.d_zy/2 + obj.restriction_y) ...
                        + obj.tail_multiplier * obj.input.footstep_plan.tail_y ...
                        - state.w_bar(2,1) / obj.input.scheme_parameters.eta ^ 2; 
            obj.y_u_M = obj.centerline_multiplier * (obj.input.footstep_plan.zmp_centerline_y + obj.input.scheme_parameters.d_zy/2 - obj.restriction_y) ...
                        + obj.tail_multiplier * obj.input.footstep_plan.tail_y ...
                        - state.w_bar(2,1) / obj.input.scheme_parameters.eta ^ 2;  
                    
            R = rotz(rad2deg(state.sf_pos(3,1))); R = R(1:2, 1:2);
            rotated_coord = R * [obj.x_u_m, obj.x_u_M; obj.y_u_m, obj.y_u_M];
            obj.x_u_m = rotated_coord(1,1);
            obj.x_u_M = rotated_coord(1,2);
            obj.y_u_m = rotated_coord(2,1);
            obj.y_u_M = rotated_coord(2,2);
            
            obj.feasibility_region(:, 1) = [obj.x_u_m; obj.x_u_M; obj.y_u_m; obj.y_u_M; 0; 0; 0; 0];
            
        end
        
        function result = getFeasibilityRegion(obj)
            
            result = obj.feasibility_region;
            
        end
        
        
    end
    
    properties (Access = private)
        
        input;
        x_u_m;
        x_u_M;
        y_u_m;
        y_u_M;
        x_u;
        y_u;
        feasibility_region;
        centerline_multiplier;
        tail_multiplier;
        prova;
        H;
        f;
        A_ineq;
        b_ineq;
        A_eq;
        b_eq;
        zmp_tracking_weight;
        P_matrix;
        restriction_x;
        restriction_y;
        restriction_builder;
        index;
        options;
        
    end
    
end