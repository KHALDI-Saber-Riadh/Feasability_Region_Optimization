%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Title: Feasability region Optimization
% Author: Saber-Riadh KHALDI, forked from: Robust Gait Generation Framework
% Author: Filippo M. Smaldone

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc ; clf ; close all; clear all;
global delta h Tc Tp TimeStep logs state plotter
delta = 0.01;
h = 0.78;
Tc = 0.7; 
Tp = 2.5;
TimeStep = 1.2;
Ltot = 0.225; % total lenght of the sole 
L_dx = 0.039; % minimum distance from the heel edges
l_bound = [0; 0.3 * Ltot];
%%%% Visualize the whole to know what type of optim to use 
% N = 10;
% l = linspace(l_bound(1), l_bound(2),N);
% L_bound = zeros(2,N);
% d_zmp = zeros(N,N);
% X = [];
% Y = [];
% tic
% for i = 1:N
%     L_bound(:,i) = [L_dx; Ltot - L_dx - l(i)];
%     L = linspace(L_bound(1, i), L_bound(2, i), N);
%     for j = 1:N
%         d_zmp(i,j) = FeasibilityRegion(L(j),l(i));
%         X = [X; l(i)];
%         Y = [Y; L(j)];
%     end
% end
% toc
% figure;
% Z = meshgrid(d_zmp);
% surf(X, Y, Z);
% xlabel('l (toe)');
% ylabel('L (Heel)');
% zlabel('d_zmp');
% title('Plot de (l, L, d_zmp)');
% title('Plot de (l, L, d_{zmp})');
%%%%%

options = optimoptions('ga', 'Display', 'off',...
                       'PopulationSize', 10, ...
                       'MaxGenerations', 15, ...
                       'MaxStallGenerations', 5);
Cost = @(x) FeasibilityRegion(x(1), x(2));
rng default  % For reproducibility
nvars = 2;
A = [0 -1;       0  1; -1 0;         1 1];
B = [ - 0.1 * Ltot  ; 0.3 * Ltot; L_dx; Ltot - L_dx]; 
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon = [];
x_opt = ga(Cost, nvars, A, B, Aeq, beq, lb, ub, nonlcon, options);
% plot the logs
% figure(2)
% for t_k = 0.1:0.1:8.0
%     plotter.plotLogsAtTimeK(logs, state, floor(t_k / delta));
% end
