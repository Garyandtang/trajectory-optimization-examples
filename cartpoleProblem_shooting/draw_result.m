%%%
% this script is to draw the result of 4 different shooting method for cart
% pole swing up problem.
%%%
% clear; clc;

load("cart_pole_result_7_1.mat")

nTrajs_list =  baseNTrajPts * (1 : size(errorResult.firstEuler, 2));

figure (55)
plot(nTrajs_list, errorResult.firstEuler);
hold on
plot(nTrajs_list, errorResult.secondEuler);
plot(nTrajs_list, errorResult.firstRk4);
plot(nTrajs_list, errorResult.secondRk4);

legend("firstEuler", "secondEuler","firstRk4","secondRk4",'Interpreter','latex','FontSize',10)