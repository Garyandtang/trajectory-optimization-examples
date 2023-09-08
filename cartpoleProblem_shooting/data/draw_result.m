%%%
% this script is to draw the result of 4 different shooting method for cart
% pole swing up problem.
%%%
clear; clc;

load("cart_pole_result_7_1.mat")

nTrajs_list =  baseNTrajPts * (1 : size(errorResult.firstEuler, 2));

figure (55)

plot(nTrajs_list, errorResult.firstEuler(1 : size(errorResult.firstEuler, 2)),":", Marker="o",LineWidth=1.5);
hold on
plot(nTrajs_list, errorResult.secondEuler(1 : size(errorResult.firstEuler, 2)), "--", Marker="o",LineWidth=1.5);
plot(nTrajs_list, errorResult.firstRk4(1 : size(errorResult.firstEuler, 2)), "-", Marker="o",LineWidth=1.5);
plot(nTrajs_list, errorResult.secondRk4(1 : size(errorResult.firstEuler, 2)), "-.", Marker="o",LineWidth=1.5);
grid on
legend("1st-Euler", "2nd-Euler","1st-RK4","2nd-RK4",'Interpreter','latex','FontSize',10)
ylabel('Total Transcription Error $\eta_{total}$','Interpreter','latex','FontSize',12)
xlabel('$N$','Interpreter','latex','FontSize',12)


figure (66)

plot(nTrajs_list, timeResult.firstEuler(1 : size(errorResult.firstEuler, 2)), ":", Marker="o",LineWidth=1.5);
hold on
plot(nTrajs_list, timeResult.secondEuler(1 : size(errorResult.firstEuler, 2)), "--", Marker="o",LineWidth=1.5);
plot(nTrajs_list, timeResult.firstRk4(1 : size(errorResult.firstEuler, 2)), "-", Marker="o",LineWidth=1.5);
plot(nTrajs_list, timeResult.secondRk4(1 : size(errorResult.firstEuler, 2)), "-.", Marker="o",LineWidth=1.5);
grid on
legend("1st-Euler", "2nd-Euler","1st-RK4","2nd-RK4",'Interpreter','latex','FontSize',10)
ylabel('NLP Solver Time $(s)$','Interpreter','latex','FontSize',12)
xlabel('$N$','Interpreter','latex','FontSize',12)


startIdx = 10;
nTrajs_list =  baseNTrajPts * (startIdx : size(errorResult.firstEuler, 2));

figure (88)
grid on
plot(nTrajs_list, errorResult.firstEuler(startIdx : size(errorResult.firstEuler, 2)), Marker="o",LineWidth=1.5);
hold on
grid on
plot(nTrajs_list, errorResult.secondEuler(startIdx : size(errorResult.firstEuler, 2)), Marker="o",LineWidth=1.5);
plot(nTrajs_list, errorResult.firstRk4(startIdx : size(errorResult.firstEuler, 2)), Marker="o",LineWidth=1.5);
plot(nTrajs_list, errorResult.secondRk4(startIdx : size(errorResult.firstEuler, 2)), Marker="o",LineWidth=1.5);
ylim([0.3 1])
