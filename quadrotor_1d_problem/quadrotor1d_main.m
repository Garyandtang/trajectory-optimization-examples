% block_move_main.m
%
% This script is used to demo the result of block move example
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Continous time trajectory formuation
%
% model:
%   see: blockMoveModel.m
%
% problem:
%       min_{u(t),q(t),dq(t)}  \int_0^1 u(t)^2 dt
%           s.t.                ddq(t) = u(t)
%                               q(0) = 0,   dq(1)=0,\\
%                               q(1) = 1,   dq(1)= 0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;

%%% Setup casadi solver 
addpath(genpath("D:\software\casadi-windows-matlabR2016a-v3.5.5"))
import casadi.*
addpath("..\models\")
addpath("..\")
% config 
config.m = 0.1;
config.g = 9.81;
%%% setup problem
model = quadrotor1d(config);
problem.model = model;

problem.cost.type = "LQR";

problem.constraints.boundarys.initState = [0;0];
problem.constraints.boundarys.finalState = [1;0];
problem.constraints.boundarys.finalTime = 1;
problem.constraints.bounds.control.lower = -inf*ones(1,1);
problem.constraints.bounds.control.upper = inf*ones(1,1);
problem.constraints.bounds.state.lower = -inf*ones(2,1);
problem.constraints.bounds.state.upper = inf*ones(2,1);
problem.grid.nTrajPts = 20;

% flag
config.flag.animationOn = false;
objApproximation = 1; % use F(t+1) - F(t) for approximation 


% first euler
config.method.dynamics = "first_order_euler";
euler1Soln = directTranscriptionMethod(problem, objApproximation, config);

% second euler
config.method.dynamics = "second_order_euler";
euler2Soln = directTranscriptionMethod(problem, objApproximation, config);

% first rk4
config.method.dynamics = "first_order_rk4";
rk1Soln = directTranscriptionMethod(problem, objApproximation, config);

% second rk4
config.method.dynamics = "second_order_rk4";
rk2Soln = directTranscriptionMethod(problem, objApproximation, config);









