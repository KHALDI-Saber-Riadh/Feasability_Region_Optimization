%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Title: Feasability region Optimization
% Author: Saber-Riadh KHALDI, forked from: Robust Gait Generation Framework
% Author: Filippo M. Smaldone

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc ; clf ; close all; clear all;
%% input data for the scheme
tic
input = struct;
L = 0.1411;
l = 0.0368;
% parameters of the scheme
input.scheme_parameters = struct;
input.scheme_parameters.delta = 0.01; % sampling time
input.scheme_parameters.h = 0.78; % CoM height for LIP model
input.scheme_parameters.g = 9.81; % gravity acceleration constant
input.scheme_parameters.eta = sqrt(input.scheme_parameters.g / input.scheme_parameters.h);
input.scheme_parameters.T_c = 0.7; % prediction horizon
input.scheme_parameters.T_p = 2.5; % preview horizon
input.scheme_parameters.C = floor(input.scheme_parameters.T_c / input.scheme_parameters.delta);
input.scheme_parameters.P = floor(input.scheme_parameters.T_p / input.scheme_parameters.delta);
input.scheme_parameters.M = 2; % optimized footstep
input.scheme_parameters.F = 4; % available footstep from the plan at each time
input.scheme_parameters.midrange = [0.15; 0.15]; % (w_mx, w_my) in m/s^2
input.scheme_parameters.dist_range = [0.15; 0.15];  % (Deltaw_mx, Deltaw_my) in m/s^2
input.scheme_parameters.dist_toe = L + 0.5 * l; % m distance between footpose and toepose
input.scheme_parameters.alpha = 0.25;
input.scheme_parameters.mi_max = 0.000025;
input.scheme_parameters.gamma_max = 0.00099;
input.scheme_parameters.epsilon = 0.0005;
input.scheme_parameters.d_zxf = L; % support polygon square width
input.scheme_parameters.d_zxb = 0.224 - L - l;
input.scheme_parameters.d_zt = l;
input.scheme_parameters.d_zy = 0.13;
input.scheme_parameters.d_ax = 0.5; % kinematic admissible region x dimension
input.scheme_parameters.d_ay = 0.25; % kinematic admissible region y dimension
input.scheme_parameters.ell_x = 0.0; % kinematic admissible region y displacement
input.scheme_parameters.ell_y = 0.2; % kinematic admissible region y displacement
input.scheme_parameters.ell = 0.2; % kinematic admissible region y displacement
input.scheme_parameters.d_ax_subsequent = 0.5; % kinematic admissible region x dimension
input.scheme_parameters.d_ay_subsequent = 0.25; % kinematic admissible region y dimension
input.scheme_parameters.ell_subsequent = 0.2; % kinematic admissible region y displacement
input.scheme_parameters.ell_x_subsequent = 0.0; % kinematic admissible region y displacement
input.scheme_parameters.ell_y_subsequent = 0.2; % kinematic admissible region y displacement
input.scheme_parameters.footstep_weight_in_cost_function = 10; % 10 %1000000;
input.scheme_parameters.zmp_track_in_cost_function = 0.01;
input.scheme_parameters.v_max = 3;
input.scheme_parameters.theta_max = 50 * pi /180;
input.scheme_parameters.omega = (2 * pi / input.scheme_parameters.T_p) * ones(input.scheme_parameters.P, 1);


%% handling non-convex constraints
% we begin with the standard convex KAR, plus two subregions that enable
% leg crossing maneuvers
input.kar = struct;
input.kar.number_of_subregions = 3;
input.kar.subregion_parameters = [input.scheme_parameters.d_ax, input.scheme_parameters.ell_x, input.scheme_parameters.d_ay, input.scheme_parameters.ell_y; ...
                                  0.15, 0.175, 0.15, 0.00; ...
                                  0.15, -0.175, 0.15, 0.00];


%% footstep plan
TimeStep = 1;
number_of_virtual_steps = 1;
input.footstep_plan = struct;
input.footstep_plan.Tp = 20;
input.footstep_plan.delta = 0.1;
P = floor(input.footstep_plan.Tp / input.footstep_plan.delta);

velocity = [0.1, 0., 0.];
input.footstep_plan.input_velocity = [velocity(1) * ones(P,1), velocity(2) * ones(P,1), velocity(3) * ones(P,1)]; % vector of input velocity over the previeuw time Vx, Vy, omega

input.footstep_plan.total_step_number = 18;
input.footstep_plan.positions = zeros(input.footstep_plan.total_step_number + number_of_virtual_steps ,3);
input.footstep_plan.timings = zeros(P ,1);
input.footstep_plan.running_steps = zeros(input.footstep_plan.total_step_number + number_of_virtual_steps ,1);
input.footstep_plan.ds_duration = 0.3 * TimeStep; % it is convenient to set a fixed duration for the double support
                                       % this can still be modified by the Step Timing Adaptation module
input.footstep_plan.dds_duration = 0.2 * TimeStep;                                       
input.footstep_plan.ds_samples = floor(input.footstep_plan.ds_duration / input.scheme_parameters.delta);
input.footstep_plan.dds_samples = floor(input.footstep_plan.dds_duration / input.scheme_parameters.delta);
input.footstep_plan.starting_sf = "right";
input.footstep_plan.tail_x = zeros(input.scheme_parameters.P - input.scheme_parameters.C, 1);
input.footstep_plan.tail_y = zeros(input.scheme_parameters.P - input.scheme_parameters.C, 1);
input.footstep_plan.zmp_centerline_x = zeros(input.scheme_parameters.C, 1);
input.footstep_plan.zmp_centerline_y = zeros(input.scheme_parameters.C, 1);
input.footstep_plan.mapping_buffer = zeros(2 * input.scheme_parameters.P, input.scheme_parameters.M + 1);
input.footstep_plan.omega = (2 * pi / input.scheme_parameters.T_p) * ones(input.scheme_parameters.P, 1);

input.sim_time = 8;

for i = 1 : P
   input.footstep_plan.timings(i, 1) = TimeStep *  (i-1);
end
input.footstep_plan.stepstimings = input.footstep_plan.timings(2:end) - input.footstep_plan.timings(1:end - 1);
input.footstep_plan.sf_pos = [0; -0.1; 0]; % position of the current support foot

input.footstep_plan.positions = [];
FSG = FootStepGenerator(input);
FSG.compute_plan();
Timing = FSG.GetFootStepTiming();
input.footstep_plan.positions = [[0, -0.1, 0]; FSG.GetFootStepPlan()];

% trick: model the initial double support as a 
% square centered between the feet
%input.footstep_plan.positions(1, 2) = 0;

% print the footstep plan
disp('input.footstep_plan.positions - (x,y,theta) [m, rad]')
disp(input.footstep_plan.positions)
% disp('input.footstep_plan.timings [s]')
% disp(input.footstep_plan.timings)



%% simulation parameters
simulation_parameters = struct;
simulation_parameters.delta = input.scheme_parameters.delta; %TODO: enable different operating frequencies
simulation_parameters.sim_time = 8;
simulation_parameters.sim_iter = 1;
simulation_parameters.sim_type = 'basic_test'; %leg_crossing'; % 'basic_test', 'leg_crossing', 'obstacle'
simulation_parameters.obstacle_number = 1;
simulation_parameters.obstacles = [0.6 -0.09 0.05 0.05]; % x,y, size_x, size_y
%% state data 
state = struct;
state.x = zeros(3,1);
state.y = zeros(3,1);
state.w_bar = zeros(2,1);
state.sf_pos = input.footstep_plan.sf_pos; % position of the current support foot
state.sf_pos_ss  = input.footstep_plan.sf_pos';
state.next_sf_pos = zeros(3,1);
state.next_sf_pos_ss = zeros(3,1);
state.current_sf = input.footstep_plan.starting_sf;
state.feasibility_region = [0; 0; 0; 0; 0; 0; 0; 0];
state.base_orient = eye(3);
state.footstep_counter = 1; % to query data from the plan
state.footstep_counter_rm = state.footstep_counter + 1;
state.footstep_counter_sm = state.footstep_counter;
state.step_time_iter = 1;
state.world_time_iter = 1;
state.sim_iter = 1;


%% log data
logs = struct;
logs.x_store = zeros(3, floor(simulation_parameters.sim_time/simulation_parameters.delta)); % x_c, x_dot_c, x_z
logs.y_store = zeros(3, floor(simulation_parameters.sim_time/simulation_parameters.delta)); % y_c, y_dot_c, y_z
logs.w_bar = zeros(2, floor(simulation_parameters.sim_time/simulation_parameters.delta)); % w_bar_x, w_bar_y
logs.w = zeros(2, floor(simulation_parameters.sim_time/simulation_parameters.delta)); % w_x, w_y
logs.actual_footsteps = zeros(3, input.footstep_plan.total_step_number); % x, y, theta
logs.feasibility_region = zeros(8, floor(simulation_parameters.sim_time/simulation_parameters.delta));
logs.feasibility_region_2 = zeros(8, floor(simulation_parameters.sim_time/simulation_parameters.delta));
logs.feasibility_region_3 = zeros(8, floor(simulation_parameters.sim_time/simulation_parameters.delta));
logs.feasibility_region_full_logs = zeros(8, floor(simulation_parameters.sim_time/simulation_parameters.delta), 6);


%% initialize plotter
plotter = Plotter(logs, input, simulation_parameters);

plotter.plotPlan();
%% CONTROLLER & SIMULATION START OPERATING HERE:
%
%
%
%
%% initialize the Robust Gait Generation Framework
wpg = RobustGaitGenerationScheme(input, state, simulation_parameters); % walking pattern generator


%% request back-up maneuver to initialize the controller (simulated)
% On a physical or multi-body simulated robot, this maneuver would consist in
% moving the CoM towards a feasible initialization for the MPC scheme:
% to do so, a simple polynomial interpolation from the current to the
% target CoM position is sufficient (to be tracked, e.g., via inverse kinematics).
% Here, we simply set the state to a proper initial value.
init_state_proposal = wpg.proposeFeasibleInitialState(state);
state.x(1,1) = init_state_proposal(1,1);
state.y(1,1) = init_state_proposal(3,1);
state.x(3,1) = state.sf_pos(1,1);
state.y(3,1) = state.sf_pos(2,1);


%% simulation cycle
for sim_iter = 1 : floor(simulation_parameters.sim_time / simulation_parameters.delta)

    % update iteration
    simulation_parameters.sim_iter = sim_iter;
 
    % solve step of gait generation algorithm
    state = wpg.update(state);
    state.sim_iter = sim_iter;

    % store logs
    logs.x_store(:, sim_iter) = state.x;
    logs.y_store(:, sim_iter) = state.y;
    logs.w_bar(:, simulation_parameters.sim_iter) = wpg.getDisturbance();
    logs.actual_footsteps(:, state.footstep_counter) = state.sf_pos;
    logs.feasibility_region(:, sim_iter) = state.feasibility_region(:, 1);
    
    for k = 1 : size(state.feasibility_region, 2)
        logs.feasibility_region_full_logs(:, sim_iter, k) = state.feasibility_region(:, k);  
    end
    

    if sim_iter == 3 && false
        plotter.plotLogs(logs, state);
    end
    if mod(sim_iter, 50) == 0 && false
        plotter.plotLogsAtTimeK(logs, state, sim_iter);
    end
    
    if state.footstep_counter >= 3
            1;    
    end  
    
    str = num2str( sim_iter / floor(simulation_parameters.sim_time / simulation_parameters.delta) * 100 );
    str = strcat('simulation is at', " ", str, '%');
    disp(str)
    
end
%
%
%
%
%% CONTROLLER & SIMULATION STOP OPERATING HERE
toc
% plot the logs
for t_k = 0.1:0.1:8.0
    plotter.plotLogsAtTimeK(logs, state, floor(t_k / input.scheme_parameters.delta));
end

zmp_x = logs.x_store(3,:);
zmp_y = logs.y_store(3,:);

com_x = logs.x_store(1,:);
com_y = logs.x_store(1,:);

diff_zmp_x = (abs(zmp_x(2:sim_iter)' - zmp_x(1:sim_iter-1)'));
diff_zmp_y = (abs(zmp_y(2:sim_iter)' - zmp_y(1:sim_iter-1)'));
d_zmp = 0;
for i = 1:sim_iter - 1
    d_zmp = d_zmp + norm([diff_zmp_x(i); diff_zmp_y(i)]);
end
