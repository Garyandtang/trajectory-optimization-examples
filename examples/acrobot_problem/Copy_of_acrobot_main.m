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
import casadi.*
addpath('utils');


% setup dynamics config 
config.m1 = 1;
config.m2 = 1;
config.g = 9.81;
config.l1 = 0.5;
config.l2 = 0.5;

%%% setup problem
problem = acrobot_swing_up(config);

% flag
config.flag.animationOn = true;
config.method.objAppro = "trapzoid_explict";

% other config
config.method.objAppro = "trapzoid_explict";
baseNTrajPts = 10;

% setup result 
timeResult.firstEuler = [];
timeResult.firstRk4 = [];
timeResult.secondEuler = [];
timeResult.secondRk4 = [];

errorResult.firstEuler = [];
errorResult.firstRk4 = [];
errorResult.secondEuler = [];
errorResult.secondRk4 = [];

% get true solution with large nTrajPts
problem.grid.nTrajPts = 500;
config.method.dynamics = "first_order_rk4";
trueSoln = directTranscriptionMethod(problem, config);

problem.trueSoln.qSoln = trueSoln.qSoln;
problem.trueSoln.dqSoln = trueSoln.dqSoln;
problem.trueSoln.ddqSoln = trueSoln.ddqSoln;
problem.trueSoln.uSoln = trueSoln.uSoln;
problem.trueSoln.tSoln = trueSoln.tSoln;

for i = 1 : 15
    %%% solve the problem
    problem.grid.nTrajPts = baseNTrajPts * i;

    % first euler
    config.method.dynamics = "first_order_euler";
    euler1Soln = directTranscriptionMethod(problem, config);

    % second euler
    config.method.dynamics = "second_order_euler";
    euler2Soln = directTranscriptionMethod(problem, config);
    
    % first rk4
    config.method.dynamics = "first_order_rk4";
    rk1Soln = directTranscriptionMethod(problem, config);
    
    % second rk4
    config.method.dynamics = "second_order_rk4";
    rk2Soln = directTranscriptionMethod(problem, config);

    %%% record the data
    % record solver time used
    timeResult.firstEuler = [timeResult.firstEuler, euler1Soln.solverTime];
    timeResult.secondEuler = [timeResult.secondEuler, euler2Soln.solverTime];
    timeResult.firstRk4 = [timeResult.firstRk4, rk1Soln.solverTime];
    timeResult.secondRk4 = [timeResult.secondRk4, rk2Soln.solverTime];

    % record error 
    errorResult.firstEuler = [errorResult.firstEuler, sum(sum(euler1Soln.info.sysDymError))];
    errorResult.secondEuler = [errorResult.secondEuler, sum(sum(euler2Soln.info.sysDymError))];
    errorResult.firstRk4 = [errorResult.firstRk4, sum(sum(rk1Soln.info.sysDymError))];
    errorResult.secondRk4 = [errorResult.secondRk4, sum(sum(rk2Soln.info.sysDymError))]
    
end

save("data\cart_pole_result.mat")





