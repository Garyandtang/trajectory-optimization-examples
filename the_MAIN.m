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
a =1
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
legend('Kelly', 'Our','Interpreter','latex')
hold off
title('System dynamic error','Interpreter','latex')
ylabel('$\varepsilon_{sd}(t)$','Interpreter','latex','FontSize',16)
xlabel('t','Interpreter','latex','FontSize',16)

% plot system dynamic error of each time interval
figure (106);clf;
plot(idx,genDynErrorIntervalSum(1,:),'o')
hold on;
plot(idx, conDynErrorIntervalSum(1,:),'o');
legend('Kelly', 'Our','Interpreter','latex')
hold off
title('System dynamic error in each time interval','Interpreter','latex','FontSize',16)
ylabel('$\eta_{sd}(t)$','Interpreter','latex','FontSize',16)
xlabel('t','Interpreter','latex','FontSize',16)

% plot system dynamic error of each time interval
figure (107);clf;
genObjCst = generalSoln.interp.objCst(t);
conObjCst = consistentSoln.interp.objCst(t);
plot(t, genObjCst(1,:));
hold on
plot(t, conObjCst(1,:));
legend('Kelly', 'Our','Interpreter','latex')
title('Objective approximation error','Interpreter','latex')
ylabel('$\varepsilon_{obj}(t)$','Interpreter','latex','FontSize',16)
xlabel('t','Interpreter','latex','FontSize',16)

return
subplot(2,2,3);
plot(t,generalCC(2,:))
hold on;
plot(t, consistentCC(2,:));
legend('general', 'Consistent','Interpreter','latex')
xlabel('time')
ylabel('d/dt pole angle')

idx = 1:length(generalSoln.info.objError);
subplot(2,2,2); hold on;
plot(idx,generalSoln.info.dynError(1,:),'bo');
plot(idx, consistentSoln.info.dynError(1,:),'ro')
hold off
legend('general', 'Consistent','Interpreter','latex')
title('State Error')
ylabel('cart position')

subplot(2,2,4); hold on;
plot(idx,generalSoln.info.dynError(2,:),'bo');
plot(idx,consistentSoln.info.dynError(2,:),'ro');
hold off
legend('general', 'Consistent','Interpreter','latex')
xlabel('segment index')
ylabel('pole angle');

figure(106); clf;
oc = generalSoln.interp.objCst(t);
subplot(1,2,1);
plot(t, oc(1,:));
title('Objective Approximation Error: ')
ylabel('d/dt cart position')

subplot(1,2,2); hold on;
plot(idx,generalSoln.info.objError(1,:),'ko');
title('State Error')
ylabel('cart position')
