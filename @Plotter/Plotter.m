classdef Plotter < handle
    
    methods (Access = public)
        
        function obj = Plotter(logs, input, simulation_parameters)
            % the constructor is just a setter
            obj.logs = logs;
            obj.sim_parameters = simulation_parameters;
            obj.input = input;
            obj.color_green = [0.4660 0.6740 0.1880];
            obj.color_blue = [0.2660 0.3740 0.8880];
            obj.dist_toe = obj.input.scheme_parameters.d_zxf + 0.5 * obj.input.scheme_parameters.d_zt;
            obj.rectangle = [obj.input.scheme_parameters.d_zxf , ...
                             obj.input.scheme_parameters.d_zxf , ...
                             -obj.input.scheme_parameters.d_zxb, ...
                             -obj.input.scheme_parameters.d_zxb, ...
                             obj.input.scheme_parameters.d_zxf ; ...
                             obj.input.scheme_parameters.d_zy / 2, ...
                             -obj.input.scheme_parameters.d_zy / 2, ...
                             -obj.input.scheme_parameters.d_zy / 2, ...
                             obj.input.scheme_parameters.d_zy / 2, ...
                             obj.input.scheme_parameters.d_zy / 2];
           obj.rectangle_toe = [obj.input.scheme_parameters.d_zt/2 , ...
                                 obj.input.scheme_parameters.d_zt/2 , ...
                                 -obj.input.scheme_parameters.d_zt/2, ...
                                 -obj.input.scheme_parameters.d_zt/2, ...
                                 obj.input.scheme_parameters.d_zt/2 ; ...
                                 obj.input.scheme_parameters.d_zy / 2, ...
                                 -obj.input.scheme_parameters.d_zy / 2, ...
                                 -obj.input.scheme_parameters.d_zy / 2, ...
                                 obj.input.scheme_parameters.d_zy / 2, ...
                                 obj.input.scheme_parameters.d_zy / 2];
            obj.initial_ds = [obj.input.scheme_parameters.d_zxf , ...
                             obj.input.scheme_parameters.d_zxf , ...
                             -obj.input.scheme_parameters.d_zxb , ...
                             -obj.input.scheme_parameters.d_zxb , ...
                             obj.input.scheme_parameters.d_zxf ; ...
                             0.1+obj.input.scheme_parameters.d_zy / 2, ...
                             -0.1-obj.input.scheme_parameters.d_zy / 2, ...
                             -0.1-obj.input.scheme_parameters.d_zy / 2, ...
                             0.1+obj.input.scheme_parameters.d_zy / 2, ...
                             0.1+obj.input.scheme_parameters.d_zy / 2];
        end
        
        function obj = plotLogs(obj, logs, state)
            
            obj.counter = obj.counter + 1;
            obj.logs = logs;
            obj.figure_handle = figure(1);
            clf;
            hold on;
            grid on;
            com = plot(logs.x_store(1, 1:state.sim_iter - 1), logs.y_store(1, 1:state.sim_iter - 1), 'r', 'Linewidth', 2);
            zmp = plot(logs.x_store(3, 1:state.sim_iter - 1), logs.y_store(3, 1:state.sim_iter - 1), 'b', 'Linewidth', 2);
            
            c = plot(obj.initial_ds(1,:), obj.initial_ds(2,:), 'm', 'Linewidth',2, 'Handlevisibility', 'off');
            
            if state.footstep_counter >= 1
                for i = 2 : state.footstep_counter
                    c = plot(obj.rectangle(1,:) + logs.actual_footsteps(1,i), obj.rectangle(2,:) + logs.actual_footsteps(2,i), 'm','Linewidth',2, 'Handlevisibility', 'off');
                    c_toe = plot(obj.rectangle_toe(1,:) + logs.actual_footsteps(1,i) + obj.dist_toe, obj.rectangle_toe(2,:) + logs.actual_footsteps(2,i), 'm','Linewidth',2, 'Handlevisibility', 'off');
                end
            end
  
            rectangle_feasibility = [state.feasibility_region(2,1), state.feasibility_region(2,1), state.feasibility_region(1,1), state.feasibility_region(1,1), state.feasibility_region(2,1); ...
                                     state.feasibility_region(4,1), state.feasibility_region(3,1), state.feasibility_region(3,1), state.feasibility_region(4,1), state.feasibility_region(4,1)];
                                 
            rectangle_feasibility_estimate = [state.feasibility_region(6,1), state.feasibility_region(6,1), state.feasibility_region(5,1), state.feasibility_region(5,1), state.feasibility_region(6,1); ...
                                     state.feasibility_region(8,1), state.feasibility_region(7,1), state.feasibility_region(7,1), state.feasibility_region(8,1), state.feasibility_region(8,1)];
                                 
            feas_region_patch = patch(rectangle_feasibility(1,:), rectangle_feasibility(2,:), obj.color_green);
            feas_region_patch.FaceAlpha = 0.3;
            feas_region_patch.EdgeColor = [0.4660 0.6740 0.1880];
            
            rectangle_feasibility = [state.feasibility_region(2,2), state.feasibility_region(2,2), state.feasibility_region(1,2), state.feasibility_region(1,2), state.feasibility_region(2,2); ...
                                     state.feasibility_region(4,2), state.feasibility_region(3,2), state.feasibility_region(3,2), state.feasibility_region(4,2), state.feasibility_region(4,2)];
            
            feas_region_patch = patch(rectangle_feasibility(1,:), rectangle_feasibility(2,:), obj.color_green);
            feas_region_patch.FaceAlpha = 0.3;
            feas_region_patch.EdgeColor = [0.2660 0.3740 0.8880];   
            
            rectangle_feasibility = [state.feasibility_region(2,3), state.feasibility_region(2,3), state.feasibility_region(1,3), state.feasibility_region(1,3), state.feasibility_region(2,3); ...
                                     state.feasibility_region(4,3), state.feasibility_region(3,3), state.feasibility_region(3,3), state.feasibility_region(4,3), state.feasibility_region(4,3)];
            
            feas_region_patch = patch(rectangle_feasibility(1,:), rectangle_feasibility(2,:), obj.color_green);
            feas_region_patch.FaceAlpha = 0.3;
            feas_region_patch.EdgeColor = [0.2660 0.3740 0.8880];     
            
            if strcmp(obj.sim_parameters.sim_type, 'obstacle')
                
                for i = 1 : obj.sim_parameters.obstacle_number
                    
                    x_m = obj.sim_parameters.obstacles(i, 1) - obj.sim_parameters.obstacles(i, 3) / 2;
                    x_M = obj.sim_parameters.obstacles(i, 1) + obj.sim_parameters.obstacles(i, 3) / 2;
                    y_m = obj.sim_parameters.obstacles(i, 2) - obj.sim_parameters.obstacles(i, 4) / 2; 
                    y_M = obj.sim_parameters.obstacles(i, 2) + obj.sim_parameters.obstacles(i, 4) / 2; 
                    object = [x_M, x_M, x_m, x_m, x_M; y_M, y_m, y_m, y_M, y_M];
                    object_patch = patch(object(1,:), object(2,:), 'k'); 
                    object_ptch.FaceAlpha = 0.3;
                    
                end
                
            end
            
%             feas_region_patch_estimate = patch(rectangle_feasibility_estimate(1,:), rectangle_feasibility_estimate(2,:), obj.color_estimate);
%             feas_region_patch_estimate.FaceAlpha = 0.1;
%             feas_region_patch_estimate.EdgeColor = [0.2660 0.3740 0.8880];            
            dcm = scatter(state.x(1,1) + state.x(2,1)/ obj.input.scheme_parameters.eta, ...
                        state.y(1,1) + state.y(2,1)/ obj.input.scheme_parameters.eta,...
                        'g','Linewidth',2);
            axis equal
            axis([-0.1 1.5 -0.25 0.55])   
            pbaspect([2 1 1]);
            legend([com, zmp, feas_region_patch, dcm], {'CoM', 'ZMP', 'current feas region', 'current dcm'});
            
        end

        function obj = plotLogsAtTimeK(obj, logs, state, k)
            
            obj.counter = obj.counter + 1;
            obj.logs = logs;
            obj.figure_handle = figure(1);
            clf;
            hold on;
            grid on;
            com = plot(logs.x_store(1, 1:k - 1), logs.y_store(1, 1:k - 1), 'r', 'Linewidth', 2);
            zmp = plot(logs.x_store(3, 1:k - 1), logs.y_store(3, 1:k - 1), 'b', 'Linewidth', 2);
            
            c = plot(obj.initial_ds(1,:), obj.initial_ds(2,:), 'm', 'Linewidth',2, 'Handlevisibility', 'off');
            
            
            ct = 0;
            for i = 1 : size(obj.input.footstep_plan.timings, 1)
                
                if obj.input.scheme_parameters.delta * k > obj.input.footstep_plan.timings(i, 1)
                    ct = i;    
                end
                
            end
            
            if ct > 1
                for i = 2 : ct
                    theta = rad2deg(logs.actual_footsteps(3,i));
                    R = rotz(theta);                    
                    dist_t = obj.dist_toe * [cosd(theta); sind(theta)];
                    center_toe = logs.actual_footsteps(1:2,i) + dist_t;
                    rot_rectangle = R(1:2, 1:2) * obj.rectangle;
                    rot_rectangle_toe = R(1:2, 1:2) * obj.rectangle_toe;
                    
                    c = plot(rot_rectangle(1,:) + logs.actual_footsteps(1,i), rot_rectangle(2,:) + logs.actual_footsteps(2,i), 'm','Linewidth',2, 'Handlevisibility', 'off');
                    c_toe = plot(rot_rectangle_toe(1,:) + center_toe(1), rot_rectangle_toe(2,:) + center_toe(2), 'm','Linewidth',2, 'Handlevisibility', 'off');
                    c_link = plot(logs.actual_footsteps(1,i), logs.actual_footsteps(2,i), 'ko','MarkerSize',8, 'MarkerFaceColor', 'k');
                end
            end
            
            for j = 1 : size(logs.feasibility_region_full_logs, 3)
                rectangle_feasibility = [logs.feasibility_region_full_logs(2, k, j), logs.feasibility_region_full_logs(2, k, j), logs.feasibility_region_full_logs(1, k, j), logs.feasibility_region_full_logs(1, k, j), logs.feasibility_region_full_logs(2, k, j); ...
                                         logs.feasibility_region_full_logs(4, k, j), logs.feasibility_region_full_logs(3, k, j), logs.feasibility_region_full_logs(3, k, j), logs.feasibility_region_full_logs(4, k, j), logs.feasibility_region_full_logs(4, k, j)];
                if j == 1
                    feas_region_patch = patch(rectangle_feasibility(1,:), rectangle_feasibility(2,:), obj.color_green);
                    feas_region_patch.FaceAlpha = 0.3;
                    feas_region_patch.EdgeColor = [0.2660 0.6740 0.1180];
                else
                    feas_region_patch = patch(rectangle_feasibility(1,:), rectangle_feasibility(2,:), obj.color_blue);
                    feas_region_patch.FaceAlpha = 0.23;
                    feas_region_patch.EdgeColor = [0.2660 0.6740 0.1180];                    
                end
            end  
            
            feas_region_patch = patch(rectangle_feasibility(1,:) + 5, rectangle_feasibility(2,:) + 5, obj.color_blue);
            feas_region_patch.FaceAlpha = 0.0;
            feas_region_patch.EdgeColor = [0.2660 0.6740 0.1180]; 
            
            if strcmp(obj.sim_parameters.sim_type, 'obstacle')
                
                for i = 1 : obj.sim_parameters.obstacle_number
                    
                    x_m = obj.sim_parameters.obstacles(i, 1) - obj.sim_parameters.obstacles(i, 3) / 2;
                    x_M = obj.sim_parameters.obstacles(i, 1) + obj.sim_parameters.obstacles(i, 3) / 2;
                    y_m = obj.sim_parameters.obstacles(i, 2) - obj.sim_parameters.obstacles(i, 4) / 2; 
                    y_M = obj.sim_parameters.obstacles(i, 2) + obj.sim_parameters.obstacles(i, 4) / 2; 
                    object = [x_M, x_M, x_m, x_m, x_M; y_M, y_m, y_m, y_M, y_M];
                    object_patch = patch(object(1,:), object(2,:), 'k'); 
                    object_ptch.FaceAlpha = 0.3;
                    
                end
                
            end
            
%             feas_region_patch_estimate = patch(rectangle_feasibility_estimate(1,:), rectangle_feasibility_estimate(2,:), obj.color_estimate);
%             feas_region_patch_estimate.FaceAlpha = 0.1;
%             feas_region_patch_estimate.EdgeColor = [0.2660 0.3740 0.8880];            
            dcm = scatter(logs.x_store(1, k) + logs.x_store(2, k) / obj.input.scheme_parameters.eta, ...
                          logs.y_store(1, k) + logs.y_store(2, k)/ obj.input.scheme_parameters.eta,...
                          'g','Linewidth',2);
            axis equal
            axis([-0.1 1.5 -0.3 0.5])
            l = legend([dcm, com, zmp, feas_region_patch], {'$x_u^k, y_u^k$','CoM', 'ZMP', 'feasibility regions'});
            if strcmp(obj.sim_parameters.sim_type, 'basic_test')
                axis([-0.1 1.5 -0.3 0.5]) 
            end
            if strcmp(obj.sim_parameters.sim_type, 'leg_crossing')
                axis([-0.1 1.5 -0.3 0.5]) 
                l = legend([dcm, com, zmp, feas_region_patch], {'$x_u^k, y_u^k$','CoM', 'ZMP', 'feasibility regions'}, 'Location', 'Southeast');
            end  
            if strcmp(obj.sim_parameters.sim_type, 'obstacle')
                axis([-0.1 1.5 -0.45 0.3]) 
                l = legend([dcm, com, zmp, feas_region_patch], {'$x_u^k, y_u^k$','CoM', 'ZMP', 'feasibility regions'}, 'Location', 'Southeast');                
            end            
            pbaspect([2 1 1]);
            
            xlabel('x [m]')
            ylabel('y [m]')
            
            set(l,'Fontsize',12)
            set(l,'Interpreter','latex')
            
        end
        
        function obj = plotPlan(obj)
            
            obj.plan_ = obj.input.footstep_plan.positions';
            obj.color_green = [0.4660 0.6740 0.1880];
            obj.color_blue = [0.2660 0.3740 0.8880];
            obj.dist_toe = obj.input.scheme_parameters.d_zxf + 0.5 * obj.input.scheme_parameters.d_zt;
            obj.rectangle = [obj.input.scheme_parameters.d_zxf , ...
                             obj.input.scheme_parameters.d_zxf , ...
                             -obj.input.scheme_parameters.d_zxb, ...
                             -obj.input.scheme_parameters.d_zxb, ...
                             obj.input.scheme_parameters.d_zxf ; ...
                             obj.input.scheme_parameters.d_zy / 2, ...
                             -obj.input.scheme_parameters.d_zy / 2, ...
                             -obj.input.scheme_parameters.d_zy / 2, ...
                             obj.input.scheme_parameters.d_zy / 2, ...
                             obj.input.scheme_parameters.d_zy / 2];
           obj.rectangle_toe = [obj.input.scheme_parameters.d_zt/2 , ...
                                 obj.input.scheme_parameters.d_zt/2 , ...
                                 -obj.input.scheme_parameters.d_zt/2, ...
                                 -obj.input.scheme_parameters.d_zt/2, ...
                                 obj.input.scheme_parameters.d_zt/2 ; ...
                                 obj.input.scheme_parameters.d_zy / 2, ...
                                 -obj.input.scheme_parameters.d_zy / 2, ...
                                 -obj.input.scheme_parameters.d_zy / 2, ... 
                                 obj.input.scheme_parameters.d_zy / 2, ...
                                 obj.input.scheme_parameters.d_zy / 2];
            obj.initial_ds = [obj.input.scheme_parameters.d_zxf , ...
                             obj.input.scheme_parameters.d_zxf , ...
                             -obj.input.scheme_parameters.d_zxb , ...
                             -obj.input.scheme_parameters.d_zxb , ...
                             obj.input.scheme_parameters.d_zxf ; ...
                             0.09+obj.input.scheme_parameters.d_zy / 2, ...
                             -0.09-obj.input.scheme_parameters.d_zy / 2, ...
                             -0.09-obj.input.scheme_parameters.d_zy / 2, ...
                             0.09+obj.input.scheme_parameters.d_zy / 2, ...
                             0.09+obj.input.scheme_parameters.d_zy / 2];
           obj.counter = obj.counter + 1;
%            obj.figure_handle = figure(1);
%            clf;
           hold on;
           grid on;
            
           c = plot(obj.initial_ds(1,:), obj.initial_ds(2,:), 'm', 'Linewidth',2, 'Handlevisibility', 'off');
            
           
           ct = length(obj.plan_);
                
               
            
           if ct > 1
               for i = 2 : ct
                    theta = rad2deg(obj.plan_(3,i));
                    R = rotz(theta);                    
                    dist_t = obj.dist_toe * [cosd(theta); sind(theta)];
                    center_toe = obj.plan_(1:2,i) + dist_t;
                    rot_rectangle = R(1:2, 1:2) * obj.rectangle;
                    rot_rectangle_toe = R(1:2, 1:2) * obj.rectangle_toe;
                    
                    c = plot(rot_rectangle(1,:) + obj.plan_(1,i), rot_rectangle(2,:) + obj.plan_(2,i), 'm','Linewidth',2, 'Handlevisibility', 'off');
                    c_toe = plot( rot_rectangle_toe(1,:) +  center_toe(1) , rot_rectangle_toe(2,:) + center_toe(2), 'm','Linewidth',2, 'Handlevisibility', 'off');
                    c_link = plot(obj.plan_(1,i), obj.plan_(2,i), 'ko','MarkerSize',8, 'MarkerFaceColor', 'k');
               end
           end
                        
           axis equal
           axis([-0.1 1.5 -0.3 0.5])
           pbaspect([2 1 1]);
            
           xlabel('x [m]')
           ylabel('y [m]')
                       
       end
        
    end
    
    properties (Access = private)
    
        logs;
        sim_parameters;
        input;
        plan_;
        figure_handle;
        rectangle;
        rectangle_toe;
        initial_ds;
        color_green;
        color_blue;
        counter;
        dist_toe;
        
    end
    
    
end
    
    