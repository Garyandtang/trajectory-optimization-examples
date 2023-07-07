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
addpath('utils');

% setup dynamics config 
config.m = 0.18;
config.g = 9.81;
config.Ixx = 0.0025;
config.l = 0.086;
%%% setup problem
problem = quadrotor_2d_move_to(config);

% flag
config.flag.animationOn = false;
config.method.objAppro = "trapzoid_explict";

% get true solution with large nTrajPts
problem.grid.nTrajPts = 500;
config.method.dynamics = "second_order_rk4";
trueSoln = directTranscriptionMethod(problem, config);

problem.trueSoln.qSoln = trueSoln.qSoln;
problem.trueSoln.dqSoln = trueSoln.dqSoln;
problem.trueSoln.ddqSoln = trueSoln.ddqSoln;
problem.trueSoln.uSoln = trueSoln.uSoln;
problem.trueSoln.tSoln = trueSoln.tSoln;

% first euler
problem.grid.nTrajPts = 30;
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

% record result
% setup result 
timeResult.firstEuler = [];
timeResult.firstRk4 = [];
timeResult.secondEuler = [];
timeResult.secondRk4 = [];

errorResult.firstEuler = [];
errorResult.firstRk4 = [];
errorResult.secondEuler = [];
errorResult.secondRk4 = [];
% record solver time used
timeResult.firstEuler = [timeResult.firstEuler, euler1Soln.solverTime];
timeResult.secondEuler = [timeResult.secondEuler, euler2Soln.solverTime];
timeResult.firstRk4 = [timeResult.firstRk4, rk1Soln.solverTime];
timeResult.secondRk4 = [timeResult.secondRk4, rk2Soln.solverTime]

% record error 
errorResult.firstEuler = [errorResult.firstEuler, sum(sum(euler1Soln.info.sysDymError))];
errorResult.secondEuler = [errorResult.secondEuler, sum(sum(euler2Soln.info.sysDymError))];
errorResult.firstRk4 = [errorResult.firstRk4, sum(sum(rk1Soln.info.sysDymError))];
errorResult.secondRk4 = [errorResult.secondRk4, sum(sum(rk2Soln.info.sysDymError))]
    
save("2d_quadrotor_shooting_results.mat")

if config.flag.animationOn
    soln = rk2Soln;
    trajhandle = @traj_diamond;
    draw_result(soln, trajhandle);
end






