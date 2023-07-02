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

% setup dynamics config 
config.m = 0.1;
config.g = 9.81;

%%% setup problem
problem = quadrotor_1d_move_to(config);

% flag
config.flag.animationOn = false;
config.method.objAppro = "trapzoid_explict";


% first euler
config.method.dynamics = "first_order_euler";
euler1Soln = directTranscriptionMethod(problem, config);

% second euler
config.method.dynamics = "second_order_euler";
euler2Soln = directTranscriptionMethod(problem,config);

% first rk4
config.method.dynamics = "first_order_rk4";
rk1Soln = directTranscriptionMethod(problem,  config);

% second rk4
config.method.dynamics = "second_order_rk4";
rk2Soln = directTranscriptionMethod(problem,  config);









