%the_MAIN.m
%
% This script is used to compare the result of general collocation methods
% and consistent collocation methods
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;

%%% provide configuration for dynamic and problem
% Number of points for initialization:
config.grid.nTrajPts = 25;

% Physical parameters for dynamics
m1 = 2; config.dyn.m1 = m1;   %cart mass
m2 = 0.1; config.dyn.m2 = m2;   %pendulum mass
config.dyn.g = 9.81;
config.dyn.l = 1;
g = config.dyn.g;
l = config.dyn.l;

% Compute an initial guess at a trajectory:
config.guess = computeGuess(config);

% Bounds:
config.bounds = computeBounds(config);
config.initState = zeros(4,1);
config.finalState = [0, pi, 0, 0]';

% flag
config.flag.animationOn = 0;

%%% Get the result 
% Kelly's Trapozoid
generalSoln = generalTrapezoidMethod(config);

% consistent trapezoid
consistentSoln = consistentTrapezoidMethod(config);

%%% Plot the results:
t = linspace(generalSoln.time(1),generalSoln.time(end),10000);
z = generalSoln.interp.state(t);
figure(101); clf; plotTraj(generalSoln,config);

figure(105); clf;
generalCC = generalSoln.interp.collCst(t);
consistentCC = consistentSoln.interp.collCst(t);

genDynErrorSum = zeros(1, size(generalCC,2));
conDynErrorSum = zeros(1, size(consistentCC,2));
genDynErrorIntervalSum = zeros(1, size(generalSoln.info.objError,2));
conDynErrorIntervalSum = zeros(1, size(consistentSoln.info.objError,2));
for i = 1:size(generalCC,1)
    genDynErrorSum = generalCC(i,:) + genDynErrorSum;
    conDynErrorSum = consistentCC(i,:) + conDynErrorSum;
    genDynErrorIntervalSum = generalSoln.info.dynError(i,:) + genDynErrorIntervalSum;
    conDynErrorIntervalSum = consistentSoln.info.dynError(i,:) + conDynErrorIntervalSum;
end
idx = 1:length(generalSoln.info.objError);

% plot system dynamic error process
figure(105); clf;
plot(t,genDynErrorSum(1,:))
hold on;
plot(t, conDynErrorSum(1,:));
legend("Kelly's", 'Ours','Interpreter','latex')
hold off
% title('System dynamic error','Interpreter','latex')
ylabel('$\varepsilon_{sd}(t)$','Interpreter','latex','FontSize',16)
xlabel('$t$','Interpreter','latex','FontSize',16)

% plot system dynamic error of each time interval
figure (106);clf;
plot(idx,genDynErrorIntervalSum(1,:),'o')
hold on;
plot(idx, conDynErrorIntervalSum(1,:),'o');
legend("Kelly's", 'Ours','Interpreter','latex')
hold off
% title('System dynamic error in each time interval','Interpreter','latex','FontSize',16)
ylabel('$\eta_{sd}(t)$','Interpreter','latex','FontSize',16)
xlabel('$t$','Interpreter','latex','FontSize',16)

% plot system dynamic error of each time interval
figure (107);clf;
genObjCst = generalSoln.interp.objCst(t);
conObjCst = consistentSoln.interp.objCst(t);
plot(t, genObjCst(1,:));
hold on
plot(t, conObjCst(1,:));
legend("Kelly's", 'Ours','Interpreter','latex')
% title('Objective approximation error','Interpreter','latex')
ylabel('$\varepsilon_{obj}(t)$','Interpreter','latex','FontSize',16)
xlabel('$t$','Interpreter','latex','FontSize',16)


