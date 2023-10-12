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
nState = problem.model.dim.nState;
nControl = problem.model.dim.nControl;
nConfig = problem.model.dim.nConfig;


% other config
config.method.objAppro = "trapzoid_implict";
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

% get true solution with large nTrajPts
problem.grid.nTrajPts = 500;
config.method.dynamics = "second_order_trapzoidal";
trueSoln = directTranscriptionMethod(problem, config);

problem.trueSoln.qSoln = trueSoln.qSoln;
problem.trueSoln.dqSoln = trueSoln.dqSoln;
problem.trueSoln.ddqSoln = trueSoln.ddqSoln;
problem.trueSoln.uSoln = trueSoln.uSoln;
problem.trueSoln.tSoln = trueSoln.tSoln;

for i = 1 : 50
    %%% solve the problem
    problem.grid.nTrajPts = baseNTrajPts * i;

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

% cp1Soln = trap1Soln;
% ourSoln = trap2Soln;
% 
% t = linspace(cp1Soln.tSoln(1),cp1Soln.tSoln(end),10000);
% idx = 1:length(ourSoln.info.sysDymError);
% 
% % plot system dynamic error process
% cp1SolnSysDymErrorTraj = cp1Soln.interp.sysDymError(t);
% ourDymErrorTraj = ourSoln.interp.sysDymError(t);
% cp1SolnSysDymErrorTrajSum = zeros(1, size(cp1SolnSysDymErrorTraj,2));
% oursDymErrorTrajSum = zeros(1, size(ourDymErrorTraj, 2));
% for i = 1 : nConfig
%     cp1SolnSysDymErrorTrajSum = cp1SolnSysDymErrorTrajSum + cp1SolnSysDymErrorTraj(i,:);
%     oursDymErrorTrajSum = oursDymErrorTrajSum + ourDymErrorTraj(i,:);
% end
% figure(105); clf;
% plot(t,cp1SolnSysDymErrorTrajSum(1,:))
% hold on;
% plot(t, oursDymErrorTrajSum(1,:));
% legend("CG1", "Ours",'Interpreter','latex','FontSize',10)
% ylabel('$\varepsilon_{dyn}(t)$','Interpreter','latex','FontSize',16)
% xlabel('$t$','Interpreter','latex','FontSize',16)
% hold off
% 
% % plot system dynamic error of each time interval
% cp1SolnSysDymErrorIntervalSum = zeros(1, size(cp1Soln.info.sysDymError,2));
% ourSolnSysDymErrorIntervalSum = zeros(1, size(ourSoln.info.sysDymError,2));
% for i = 1 : nConfig
%     cp1SolnSysDymErrorIntervalSum = cp1SolnSysDymErrorIntervalSum + cp1Soln.info.sysDymError(i,:);
%     ourSolnSysDymErrorIntervalSum = ourSolnSysDymErrorIntervalSum + ourSoln.info.sysDymError(i,:);
% end
% figure(106); clf;
% 
% scatter(idx,cp1SolnSysDymErrorIntervalSum(1,:),'o','filled')
% box on;
% hold on;
% scatter(idx, ourSolnSysDymErrorIntervalSum(1,:),'o');
% legend("CG1","Ours",'Interpreter','latex')
% ylabel('$\eta_{dyn,k}$','Interpreter','latex','FontSize',16)
% xlabel('$k$','Interpreter','latex','FontSize',16)
% hold off
%  sum(cp1SolnSysDymErrorIntervalSum)
%  sum(ourSolnSysDymErrorIntervalSum)

