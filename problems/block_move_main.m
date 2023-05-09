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

%%% setup problem
model = blockMoveModel();
problem.model = model;

problem.cost.type = "LQR";

problem.constraints.boundarys.initState = zeros(2,1);
problem.constraints.boundarys.finalState = [1;0];
problem.constraints.boundarys.finalTime = 1;
problem.constraints.bounds.control.lower = -inf*ones(1,1);
problem.constraints.bounds.control.upper = inf*ones(1,1);
problem.constraints.bounds.state.lower = -inf*ones(2,1);
problem.constraints.bounds.state.upper = inf*ones(2,1);
problem.grid.nTrajPts = 30;

firstOrderSoln = firstOrderTrapzoidMethod(problem)
second = consistentTrapzoidMethod(problem)