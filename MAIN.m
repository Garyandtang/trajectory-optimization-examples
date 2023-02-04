% MAIN.m
%
% This script runs trajectory optimization using either orthogonal
% collocation or direct transcription (trapazoid).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem description
%
% Cartpole swing up: 
%   State:
%       z = [4, 1] = state vector = [x;q;dx;dq]
%   
%   Control:
%       u = [1, 1] = actuation vector = F = force on cart
%
%   objective:
%       L = u^2
%   Parameters:
%       .m1 = cart mass
%       .m2 = pendulum point-mass
%       .g = gravity
%       .l = length of the pendulum
%     
%   System dynamics:
%       dz = f(z, u)
%       dz(1) = z(3)  % first order
%       dz(2) = z(4)  % first order
%       dz(3) = ddx
%       dz(4) = ddq
%       ddx = l*m2*sin(q)*dq2^2 + u + m2*g*cos(q)*sin(q)/(m1+m2*(1 - cos(q)^2))
%       ddq = -(l*m2*cos(q)*sin(q)*dq2^2 + u*cos(q)+(m1+m2)*g*sin(q))/(l*m1+l*m2*(1 - cos(q)^2))
%
%
%   Constraints:
%       -49 <= u <= 49      bounded control input 
%       -1 <= x <= 1        bounded position
%       z_init = [0,0,0,0]  init state constraint
%       z_f = [0,pi,0,0]    final state constraint
%
%   No of control input:
%       N = 15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear;


% Number of points for initialization:
config.grid.nTrajPts = 15;

% Options for nlp solver (fmincon)
config.options.nlp = optimset(...
    'Display','iter',...
    'MaxIter',250,...
    'MaxFunEvals',5e4);

% Physical parameters for dynamics
m1 = 1.0; config.dyn.m1 = m1;   %cart mass
m2 = 1.0; config.dyn.m2 = m2;   %pendulum mass
config.dyn.g = 9.81;
config.dyn.l = 1;

% Compute an initial guess at a trajectory:
config.guess = computeGuess(config);

% Bounds:
config.bounds = computeBounds(config);

% Create function handles to be called by optimization:
config.function.dynamics = @(t,z,u)( cartPoleDynamics(z,u,config.dyn) );
config.function.costIntegrand = @(t,z,u)( costFunction(t,z,u) );

%%%% Select method and compute trajectory
% traj = orthogonalCollocation(config);
traj = dirTrans_Trapz(config);

% Animation:
P.plotFunc = @(t,z)( drawCartPole(t,z,config.dyn) );
P.speed = 0.7;
P.figNum = 102;
t = linspace(traj.time(1),traj.time(end),250);
z = traj.interp.state(t);
animate(t,z,P)

% Plot the results:
figure(101); clf; plotTraj(traj,config);