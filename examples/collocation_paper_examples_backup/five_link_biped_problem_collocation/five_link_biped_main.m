clc;clear;

%%% Setup casadi solver 
import casadi.*
addpath('utils');

% setup dynamics config 
config.m1 = 3.2;
config.m2 = 6.8;
config.m3 = 10;
config.m4 = 6.8;
config.m5 = 3.2;
config.I1 = 0.93;
config.I2 = 1.0800;
config.I3 = 2.2200;
config.I4 = 1.0800;
config.I5 = 0.9300;
config.l1 = 0.4000;
config.l2 = 0.4000;
config.l3 = 0.6250;
config.l4 = 0.4000;
config.l5 = 0.4000;
config.c1 = 0.1280;
config.c2 = 0.1630;
config.c3 = 0.2000;
config.c4 = 0.1630;
config.c5 = 0.1280;
config.g = 9.8100;
config.stepLength = 0.5000;
config.stepTime = 0.7000;

%%% setup problem
problem = five_link_biped_walk(config);

% flag
config.flag.animationOn = true;
config.method.objAppro = "trapzoid_implict";

% get true solution with large nTrajPts
problem.grid.nTrajPts = 500;
config.method.dynamics = "second_order_trapzoidal";
trueSoln = directTranscriptionMethod(problem, config);

problem.trueSoln.qSoln = trueSoln.qSoln;
problem.trueSoln.dqSoln = trueSoln.dqSoln;
problem.trueSoln.ddqSoln = trueSoln.ddqSoln;
problem.trueSoln.uSoln = trueSoln.uSoln;
problem.trueSoln.tSoln = trueSoln.tSoln;
baseNTrajPts = 5;
% setup result 
timeResult.cg1 = [];
timeResult.cg2 = [];
timeResult.cg3 = [];
timeResult.ours = [];

errorResult.cg1 = [];
errorResult.cg2 = [];
errorResult.cg3 = [];
errorResult.ours = [];
for i = 1 : 12
    %%% solve the problem
    problem.grid.nTrajPts = baseNTrajPts * i ;

    % first euler
    config.method.dynamics = "first_order_trapzoidal";
    config.method.objAppro = "trapzoid_explict";
    cg1Soln = directTranscriptionMethod(problem, config);

    % second euler
    config.method.dynamics = "first_order_trapzoidal";
    config.method.objAppro = "trapzoid_implict";
    cg2Soln = directTranscriptionMethod(problem, config);
    
    % first rk4
    config.method.dynamics = "second_order_trapzoidal";
    config.method.objAppro = "trapzoid_explict";
    cg3Soln = directTranscriptionMethod(problem, config);
    
    % second rk4
    config.method.dynamics = "second_order_trapzoidal";
    config.method.objAppro = "trapzoid_implict";
    oursSoln = directTranscriptionMethod(problem, config);

    %%% record the data
    % record solver time used
    timeResult.cg1 = [timeResult.cg1, cg1Soln.solverTime];
    timeResult.cg2 = [timeResult.cg2, cg2Soln.solverTime];
    timeResult.cg3 = [timeResult.cg3, cg3Soln.solverTime];
    timeResult.ours = [timeResult.ours, oursSoln.solverTime];

    % record error 
    errorResult.cg1 = [errorResult.cg1, sum(sum(cg1Soln.info.sysDymError))];
    errorResult.cg2 = [errorResult.cg2, sum(sum(cg2Soln.info.sysDymError))];
    errorResult.cg3 = [errorResult.cg3, sum(sum(cg3Soln.info.sysDymError))];
    errorResult.ours = [errorResult.ours, sum(sum(oursSoln.info.sysDymError))]
    
end



save("data\cart_pole_result.mat")




