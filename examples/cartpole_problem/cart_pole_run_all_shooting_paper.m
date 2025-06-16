clc;clear;

%%% Setup casadi solver 
import casadi.*

%%% setup dynamics config
config.dyn.g = 9.81;
config.dyn.l = 1;
config.dyn.m1 = 1;
config.dyn.m2 = 0.1;

%%% setup problem
problem = cart_pole_swing_up(config);

% other config
config.method.objAppro = "trapzoid_implict";
baseNTrajPts = 5;

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
problem.grid.nTrajPts = 150;
config.method.dynamics = "second_order_rk4";
trueSoln = directTranscriptionMethod(problem, config);

problem.trueSoln.qSoln = trueSoln.qSoln;
problem.trueSoln.dqSoln = trueSoln.dqSoln;
problem.trueSoln.ddqSoln = trueSoln.ddqSoln;
problem.trueSoln.uSoln = trueSoln.uSoln;
problem.trueSoln.tSoln = trueSoln.tSoln;

for i = 1 : 12
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


