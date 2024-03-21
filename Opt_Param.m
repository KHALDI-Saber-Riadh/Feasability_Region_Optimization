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
N = 1;
l = linspace(l_bound(1), l_bound(2),N);
L_bound = zeros(2,N);
d_zmp = zeros(N,N);
tic
for i = 1:N
    L_bound(:,i) = [L_dx; Ltot - L_dx - l(i)];
    L = linspace(L_bound(1), L_bound(2), N);
    for j = 1:N
        d_zmp(i,j) = FeasibilityRegion(L(j),l(i));
    end
end
toc
% Cost = @(x) FeasibilityRegion(x(1), x(2));
% rng default  % For reproducibility
% nvars = 2;
% lb = [L_bound(1) , l_bound(1)];
% ub = [L_bound(2), l_bound(2)];
% x = particleswarm(Cost,nvars,lb,ub);

% plot the logs
for t_k = 0.1:0.1:8.0
    plotter.plotLogsAtTimeK(logs, state, floor(t_k / delta));
end