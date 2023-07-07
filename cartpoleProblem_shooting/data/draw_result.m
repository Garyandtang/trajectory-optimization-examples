%%%
% this script is to draw the result of 4 different shooting method for cart
% pole swing up problem.
%%%
clear; clc;

load("cart_pole_result_7_1.mat")

nTrajs_list =  baseNTrajPts * (1 : size(errorResult.firstEuler, 2));

figure (55)
plot(nTrajs_list, errorResult.firstEuler(1 : size(errorResult.firstEuler, 2)), Marker="o");
hold on
plot(nTrajs_list, errorResult.secondEuler(1 : size(errorResult.firstEuler, 2)), Marker="o");
plot(nTrajs_list, errorResult.firstRk4(1 : size(errorResult.firstEuler, 2)), Marker="o");
plot(nTrajs_list, errorResult.secondRk4(1 : size(errorResult.firstEuler, 2)), Marker="o");

legend("Euler1", "Euler2","RK1","RK2",'Interpreter','latex','FontSize',10)
ylabel('Total Transcription Error $\eta_{total}$','Interpreter','latex','FontSize',12)
xlabel('$N$','Interpreter','latex','FontSize',12)


figure (66)
plot(nTrajs_list, timeResult.firstEuler(1 : size(errorResult.firstEuler, 2)), Marker="o");
hold on
plot(nTrajs_list, timeResult.secondEuler(1 : size(errorResult.firstEuler, 2)), Marker="o");
plot(nTrajs_list, timeResult.firstRk4(1 : size(errorResult.firstEuler, 2)), Marker="o");
plot(nTrajs_list, timeResult.secondRk4(1 : size(errorResult.firstEuler, 2)), Marker="o");

legend("Euler1", "Euler2","RK1","RK2",'Interpreter','latex','FontSize',10)
ylabel('NLP Solver Time $(s)$','Interpreter','latex','FontSize',12)
xlabel('$N$','Interpreter','latex','FontSize',12)