% clc;clear;

%%% Setup casadi solver 
addpath(genpath("D:\software\casadi-windows-matlabR2016a-v3.5.5"))
import casadi.*
addpath("..\models\")
addpath("..\");
addpath("..\solvers\")

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
problem.grid.nTrajPts = 100;
baseNTrajPts = 10;
% flag
config.flag.animationOn = false;
config.method.objAppro = "explict_trapzoid";
objApproximation = 1; % use F(t+1) - F(t) for approximation 

% setup result 
timeResult.firstEuler = [];
timeResult.firstRk4 = [];
timeResult.secondEuler = [];
timeResult.secondRk4 = [];

errorResult.firstEuler = [];
errorResult.firstRk4 = [];
errorResult.secondEuler = [];
errorResult.secondRk4 = [];

% % get true solution with large nTrajPts
% problem.grid.nTrajPts = 500;
% config.method.dynamics = "second_order_euler";
% trueSoln = directTranscriptionMethod(problem, objApproximation, config);
% 
% problem.trueSoln.qSoln = trueSoln.qSoln;
% problem.trueSoln.dqSoln = trueSoln.dqSoln;
% problem.trueSoln.ddqSoln = trueSoln.ddqSoln;
% problem.trueSoln.uSoln = trueSoln.uSoln;
% problem.trueSoln.tSoln = trueSoln.tSoln;

for i = 1 : 15
    %%% solve the problem
    problem.grid.nTrajPts = baseNTrajPts * i;

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

save("cart_pole_result_7_1.mat")


return

%%% solve the problem and get solution
config.method.objAppro = "explict_trapzoid";
config.method.dynamics = "first_order_rk4";
objApproximation = 1; % use F(t+1) - F(t) for approximation 
cp1Soln = firstOrderShooting(problem,objApproximation, config);

return
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

% t = linspace(euler1Soln.tSoln(1),cp1Soln.tSoln(end),1000);
% idx = 1:length(euler1Soln.info.sysDymError);
% 
% euler1Error = zeros(1,size(euler1Soln.info.sysDymError,2));
% euler2Error = zeros(1,size(euler1Soln.info.sysDymError,2));
% rk1Error = zeros(1,size(euler1Soln.info.sysDymError,2));
% rk2Error = zeros(1,size(euler1Soln.info.sysDymError,2));
% 
% for j = 1 : nControl
%     euler1Error = euler1Error + euler1Soln.info.sysDymError(j, :);
%     euler2Error = euler2Error + euler2Soln.info.sysDymError(j, :);
%     rk1Error = rk1Error + rk1Soln.info.sysDymError(j, :);
%     rk2Error = rk2Error + rk2Soln.info.sysDymError(j, :);
% end