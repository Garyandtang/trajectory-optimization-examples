clc;clear;

%%% Setup casadi solver 
addpath(genpath("D:\software\casadi-windows-matlabR2016a-v3.5.5"))
import casadi.*
addpath("..\models\")
addpath("..\");
addpath("..\solvers\")

%%% setup dynamics config
config.dyn.g = 9.81;
config.dyn.l = 0.5;
config.dyn.m1 = 1;
config.dyn.m2 = 0.5;



%%% setup problem
problem = cart_pole_swing_up(config);

problem.constraints.boundarys.finalState = [0;pi;0;0];
problem.constraints.boundarys.finalTime = 3;
baseNTrajPts = 30;




% % get true solution with large nTrajPts
% problem.grid.nTrajPts = 1000;
% config.method.dynamics = "second_order_trapzoidal";
% config.method.objAppro = "trapzoid_implict";
% trueSoln = directTranscriptionMethod(problem, config);
% 
% problem.trueSoln.qSoln = trueSoln.qSoln;
% problem.trueSoln.dqSoln = trueSoln.dqSoln;
% problem.trueSoln.ddqSoln = trueSoln.ddqSoln;
% problem.trueSoln.uSoln = trueSoln.uSoln;
% problem.trueSoln.tSoln = trueSoln.tSoln;

problem.grid.nTrajPts = baseNTrajPts;
% CG1
config.method.dynamics = "first_order_trapzoidal";
config.method.objAppro = "trapzoid_explict";
cg1Soln = directTranscriptionMethod(problem, config);

% CG2
config.method.dynamics = "first_order_trapzoidal";
config.method.objAppro = "trapzoid_implict";
cg2Soln = directTranscriptionMethod(problem, config);

% CG3
config.method.dynamics = "second_order_trapzoidal";
config.method.objAppro = "trapzoid_explict";
cg3Soln = directTranscriptionMethod(problem, config);

% Ours
config.method.dynamics = "second_order_trapzoidal";
config.method.objAppro = "trapzoid_implict";
oursSoln = directTranscriptionMethod(problem, config);


t = linspace(oursSoln.tSoln(1),oursSoln.tSoln(end),5000);
cg1SysDynError = cg1Soln.interp.sysDymError(t);
cg2SysDynError = cg2Soln.interp.sysDymError(t);
cg3SysDynError = cg3Soln.interp.sysDymError(t);
oursSysDynError = oursSoln.interp.sysDymError(t);

% plot system dynamics error
figure(65); clf;
ele = 2;
plot(t,abs(cg1SysDynError(ele,:)),'-',LineWidth=1.5);
box on;
grid on
hold on

plot(t,abs(cg2SysDynError(ele,:)),'--',LineWidth=1.5);
plot(t,abs(cg3SysDynError(ele,:)),':',LineWidth=1.5);
plot(t,abs(oursSysDynError(ele,:)),'-.',LineWidth=1.5);
legend("CG1", 'CG2', "CG3","Ours",'Interpreter','latex','FontSize',10)
ylabel('Approximation Error $\varepsilon(t)$','Interpreter','latex','FontSize',16)
xlabel('$t$','Interpreter','latex','FontSize',16)


