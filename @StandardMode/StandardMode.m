classdef StandardMode < FeasibilityDrivenBase
    
    methods (Access = public)
        
        function obj = StandardMode(input)
            
            obj.input = input;
            obj.x_u_m = 0;
            obj.x_u_M = 0;
            obj.y_u_m = 0;
            obj.y_u_M = 0;
            obj.feasibility_region = [obj.x_u_m; obj.x_u_M; obj.y_u_m; obj.y_u_M];
            obj.centerline_multiplier = obj.input.scheme_parameters.delta * ...
                                        exp(obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta * (0 : obj.input.scheme_parameters.C - 1));
            obj.tail_multiplier = obj.input.scheme_parameters.delta * ...
                                        exp(obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta * (obj.input.scheme_parameters.C : obj.input.scheme_parameters.P - 1));            
        end
        
        function obj = solve(obj)
            
        end
        
        function is_feasible = feasibilityCheck(obj, state, input)
            obj.input = input;
            obj.x_u_m = obj.centerline_multiplier * (obj.input.footstep_plan.zmp_centerline_x - obj.input.scheme_parameters.d_z/2) ...
                        + obj.tail_multiplier * obj.input.footstep_plan.tail_x ...
                        - state.w_bar(1,1) / obj.input.scheme_parameters.eta ^ 2;
            obj.x_u_M = obj.centerline_multiplier * (obj.input.footstep_plan.zmp_centerline_x + obj.input.scheme_parameters.d_z/2) ...
                        + obj.tail_multiplier * obj.input.footstep_plan.tail_x ...
                        - state.w_bar(1,1) / obj.input.scheme_parameters.eta ^ 2; 
            obj.y_u_m = obj.centerline_multiplier * (obj.input.footstep_plan.zmp_centerline_y - obj.input.scheme_parameters.d_z/2) ...
                        + obj.tail_multiplier * obj.input.footstep_plan.tail_y ...
                        - state.w_bar(2,1) / obj.input.scheme_parameters.eta ^ 2; 
            obj.y_u_M = obj.centerline_multiplier * (obj.input.footstep_plan.zmp_centerline_y + obj.input.scheme_parameters.d_z/2) ...
                        + obj.tail_multiplier * obj.input.footstep_plan.tail_y ...
                        - state.w_bar(2,1) / obj.input.scheme_parameters.eta ^ 2;  
            obj.feasibility_region = [obj.x_u_m; obj.x_u_M; obj.y_u_m; obj.y_u_M];
            
            is_feasible = false;
            obj.x_u = state.x(1,1) + state.x(2,1) / obj.input.scheme_parameters.eta;
            obj.y_u = state.y(1,1) + state.y(2,1) / obj.input.scheme_parameters.eta; 
            if (obj.x_u >= obj.x_u_m && obj.x_u <= obj.x_u_M) && (obj.y_u >= obj.y_u_m && obj.y_u <= obj.y_u_M)
                is_feasible = true;    
            end
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
    end
    
    %properties (Access = protected)
    %    input = 1;
    %end
    
end