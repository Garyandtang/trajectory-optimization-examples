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
addpath("..")
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
problem.grid.nTrajPts = 20;


%%% solve the problem and get solution
% case 1 with obj appro
objApproximation = 1; % use F(t+1) - F(t) for approximation 
cp2Soln = firstOrderTrapzoidMethod(problem,objApproximation);
consistentSoln = consistentTrapzoidMethod(problem, objApproximation);

% case 2 without obj appro
objApproximation = 0; % use trapzoid method for obj approximation
cp1Soln = firstOrderTrapzoidMethod(problem, objApproximation);
cp3Soln = consistentTrapzoidMethod(problem, objApproximation);

%%% Plot the results:
t = linspace(consistentSoln.tSoln(1),consistentSoln.tSoln(end),1000);
idx = 1:length(consistentSoln.info.configError);
consistentQTraj = consistentSoln.interp.q(t);
optimalQTraj = consistentSoln.optimalSoln.q(t);

% plot result q trajectory and optimal q trajectory
figure(101); clf; plot(t,consistentQTraj);
hold on
plot(t, optimalQTraj);

% plot configuration error evolution
figure(65); clf;
plot(t,cp1Soln.interp.configError(t),'-',LineWidth=1.5);
box on;
grid on
hold on
plot(t,cp2Soln.interp.configError(t),'--',LineWidth=1.5);
plot(t,cp3Soln.interp.configError(t),':',LineWidth=1.5);
plot(t,consistentSoln.interp.configError(t),'-.',LineWidth=1.5);
legend("CG1", 'CG2', "CG3","Ours",'Interpreter','latex','FontSize',10)
ylabel('Approximation Error $\varepsilon(t)$','Interpreter','latex','FontSize',16)
xlabel('$t$','Interpreter','latex','FontSize',16)
% plot configuration error in each time interval
figure(66); clf; 
scatter(idx,cp1Soln.info.configError(1,:), 'o','filled');
box on;
grid on
hold on

scatter(idx,cp2Soln.info.configError(1,:), '^','filled');
scatter(idx,cp3Soln.info.configError(1,:), 'v','filled');
scatter(idx,consistentSoln.info.configError(1,:), 'square','filled');
legend("CG1", 'CG2', "CG3","Ours",'Interpreter','latex')
ylabel('Accumulated Error $\eta_k$','Interpreter','latex','FontSize',16)
xlabel('$k$','Interpreter','latex','FontSize',16)
hold off
% plot system dynamic error evolution
figure(67); clf;
cp1SysDymError = cp1Soln.interp.sysDymError(t);
plot(t,cp1SysDymError(end,:));
hold on
cp2SysDymError = cp2Soln.interp.sysDymError(t);
plot(t,cp2SysDymError(end,:));
% hold on
cp3SysDymError = cp3Soln.interp.sysDymError(t);
plot(t,cp3SysDymError(end,:));
plot(t,consistentSoln.interp.sysDymError(t));
legend("CG1", 'CG2', "CG3","Ours",'Interpreter','latex')
