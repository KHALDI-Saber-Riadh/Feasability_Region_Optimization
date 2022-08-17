classdef RobustGaitGenerationScheme

    methods (Access = public)

        function obj = RobustGaitGenerationScheme(input, state)
            % constructor
            obj.input = input;
            obj.sm_instance = StandardMode(obj.input);
            obj.rm_instance = RecoveryMode(obj.input);
            obj.dob_instance = DisturbanceObserver(obj.input, state);
            obj.restriction_builder = BuildRestrictionFunction(obj.input);
            
            obj.steps_in_horizon = zeros(obj.input.scheme_parameters.F + 2, 2);   
            obj.centerline_temp_temp_x = zeros(2 * obj.input.scheme_parameters.P, 1);
            obj.centerline_temp_temp_y = zeros(2 * obj.input.scheme_parameters.P, 1); 
            obj.centerline_temp_x = cell(obj.input.scheme_parameters.F + 1, 1);
            obj.centerline_temp_y = cell(obj.input.scheme_parameters.F + 1, 1);
            
        end

        function obj = update(obj, state)
            
             % detect obstacles
             % TODO
             
             % observer
             zmpdot = zeros(2,1);
             obj.dob_instance.update([state.x(1,1); state.x(3,1)], ...
                                     [state.y(1,1); state.y(3,1)], ...
                                     zmpdot);  
             state.w_bar = obj.dob_instance.getDisturbance();
             
             
             % build a ZMP trajectory to be used when needed                                 
             obj.steps_in_horizon(1 : obj.input.scheme_parameters.F + 2, :) = ...
                 obj.input.footstep_plan.positions(state.footstep_counter : state.footstep_counter + obj.input.scheme_parameters.F + 1, 1:2);
             
             time_counter = 1;
             for i = 1 : obj.input.scheme_parameters.F + 1
                 
                 obj.ss_samples = round( (obj.input.footstep_plan.timings(i + 1, 1) - obj.input.footstep_plan.timings(i, 1)) / obj.input.scheme_parameters.delta ) - ...
                                  obj.input.footstep_plan.ds_samples;
                 obj.centerline_temp_x{i} = [obj.steps_in_horizon(i, 1) + (obj.steps_in_horizon(i + 1, 1) - obj.steps_in_horizon(i, 1)) * ...
                                            (0 : obj.input.footstep_plan.ds_samples - 1)' / obj.input.footstep_plan.ds_samples; ...
                                            obj.steps_in_horizon(i + 1, 1) * ones(obj.ss_samples ,1)];
                 obj.centerline_temp_y{i} = [obj.steps_in_horizon(i, 2) + (obj.steps_in_horizon(i + 1, 2) - obj.steps_in_horizon(i, 2)) * ...
                                            (0 : obj.input.footstep_plan.ds_samples - 1)' / obj.input.footstep_plan.ds_samples; ...
                                            obj.steps_in_horizon(i + 1, 2) * ones(obj.ss_samples ,1)]; 
                 
                 obj.centerline_temp_temp_x(time_counter :  time_counter + obj.ss_samples + obj.input.footstep_plan.ds_samples - 1, 1) = obj.centerline_temp_x{i};                    
                 obj.centerline_temp_temp_y(time_counter :  time_counter + obj.ss_samples + obj.input.footstep_plan.ds_samples - 1, 1) = obj.centerline_temp_y{i};                  
                 time_counter = time_counter + obj.ss_samples + obj.input.footstep_plan.ds_samples;
                 
                 
             end
             
             index = state.world_time_iter - round( obj.input.footstep_plan.timings(state.footstep_counter, 1) / obj.input.scheme_parameters.delta );
             obj.input.footstep_plan.zmp_centerline_x = obj.centerline_temp_temp_x(index : index + obj.input.scheme_parameters.C - 1, 1);
             obj.input.footstep_plan.zmp_centerline_y = obj.centerline_temp_temp_y(index : index + obj.input.scheme_parameters.C - 1, 1);
             obj.input.footstep_plan.tail_x = obj.centerline_temp_temp_x(index + obj.input.scheme_parameters.C: index + obj.input.scheme_parameters.P - 1, 1);                                                                           
             obj.input.footstep_plan.tail_y = obj.centerline_temp_temp_y(index + obj.input.scheme_parameters.C: index + obj.input.scheme_parameters.P - 1, 1);                                                                                                
             
             
             % feasibility check (only if cleared by the obstacle
             % detection)
             obj.sm_instance.feasibilityCheck(state, obj.input);
               
             % maybe STA (?) 
               % recovery mode . STA
    
               
             % control cycle
             
             
             % integrate LIP
             obj.integrateModel();
             
             
        end

        function w_bar = getDisturbance(obj)

            w_bar = obj.dob_instance.getDisturbance();

        end
        
        function data = getModifiedInputStructure(obj)

            data = obj.input;

        end

    end
    
    methods (Access = private) 
       
        function obj = integrateModel(obj)
            
            
        end
        
    end

    properties (Access = private)
        
        % building blocks
        input;
        sm_instance;
        rm_instance;
        dob_instance;
        restriction_builder;
        
        % buffers
        steps_in_horizon;
        centerline_temp_x;
        centerline_temp_y;
        centerline_temp_temp_x;
        centerline_temp_temp_y;      
        ss_samples;
        
    end
    
end