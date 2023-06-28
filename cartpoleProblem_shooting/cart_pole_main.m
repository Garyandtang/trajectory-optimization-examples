clc;clear;

%%% Setup casadi solver 
addpath(genpath("D:\software\casadi-windows-matlabR2016a-v3.5.5"))
import casadi.*
addpath("..\models\")
addpath("../");
%%% setup problem
config.dyn.g = 9.81;
config.dyn.l = 1;
config.dyn.m1 = 1;
config.dyn.m2 = 0.1;
model = cartPoleModel(config);
problem.model = model;
nState = model.dim.nState;
nControl = model.dim.nControl;
nConfig = model.dim.nConfig;
problem.cost.type = "LQR";
% wn = sqrt(l/g);
problem.constraints.boundarys.initState = zeros(nState,1);
problem.constraints.boundarys.finalState = [0;pi;0;0];
problem.constraints.boundarys.finalTime = 3;
problem.constraints.bounds.control.lower = -inf*ones(1,1);
problem.constraints.bounds.control.upper = inf*ones(1,1);
problem.constraints.bounds.state.lower = -inf*ones(4,1);
problem.constraints.bounds.state.upper = inf*ones(4,1);
problem.grid.nTrajPts = 10;

% setup result 
timeResult.firstEuler = [];
timeResult.firstRk4 = [];
timeResult.secondEuler = [];
timeResult.secondRk4 = [];

errorResult.firstEuler = [];
errorResult.firstRk4 = [];
errorResult.secondEuler = [];
errorResult.secondRk4 = [];


for i = 1 : 10
end

% flag
config.flag.animationOn = false;



%%% solve the problem and get solution
config.method.objAppro = "explict_trapzoid";
config.method.dynamics = "first_order_rk4";
objApproximation = 1; % use F(t+1) - F(t) for approximation 
cp1Soln = firstOrderShooting(problem,objApproximation, config);

config.method.objAppro = "explict_trapzoid";
config.method.dynamics = "second_order_rk4";
ourSoln = firstOrderShooting(problem,objApproximation, config);



t = linspace(cp1Soln.tSoln(1),cp1Soln.tSoln(end),10000);
idx = 1:length(ourSoln.info.sysDymError);

% plot system dynamic error process
cp1SolnSysDymErrorTraj = cp1Soln.interp.sysDymError(t);
ourDymErrorTraj = ourSoln.interp.sysDymError(t);
cp1SolnSysDymErrorTrajSum = zeros(1, size(cp1SolnSysDymErrorTraj,2));
oursDymErrorTrajSum = zeros(1, size(ourDymErrorTraj, 2));
for i = 1 : nConfig
    cp1SolnSysDymErrorTrajSum = cp1SolnSysDymErrorTrajSum + cp1SolnSysDymErrorTraj(i,:);
    oursDymErrorTrajSum = oursDymErrorTrajSum + ourDymErrorTraj(i,:);
end
figure(105); clf;
plot(t,cp1SolnSysDymErrorTrajSum(1,:))
hold on;
plot(t, oursDymErrorTrajSum(1,:));
legend("CG1", "Ours",'Interpreter','latex','FontSize',10)
ylabel('$\varepsilon_{dyn}(t)$','Interpreter','latex','FontSize',16)
xlabel('$t$','Interpreter','latex','FontSize',16)
hold off

% plot system dynamic error of each time interval
cp1SolnSysDymErrorIntervalSum = zeros(1, size(cp1Soln.info.sysDymError,2));
ourSolnSysDymErrorIntervalSum = zeros(1, size(ourSoln.info.sysDymError,2));
for i = 1 : nConfig
    cp1SolnSysDymErrorIntervalSum = cp1SolnSysDymErrorIntervalSum + cp1Soln.info.sysDymError(i,:);
    ourSolnSysDymErrorIntervalSum = ourSolnSysDymErrorIntervalSum + ourSoln.info.sysDymError(i,:);
end
figure(106); clf;

scatter(idx,cp1SolnSysDymErrorIntervalSum(1,:),'o','filled')
box on;
hold on;
scatter(idx, ourSolnSysDymErrorIntervalSum(1,:),'o');
legend("CG1","Ours",'Interpreter','latex')
ylabel('$\eta_{dyn,k}$','Interpreter','latex','FontSize',16)
xlabel('$k$','Interpreter','latex','FontSize',16)
hold off

sum(cp1SolnSysDymErrorIntervalSum)
sum(ourSolnSysDymErrorIntervalSum)