%%%
% this script is to draw the result of 4 different shooting method for cart
% pole swing up problem.
%%%
clear; clc;

load("cart_pole_result.mat")
startIdx = 2;
endIdx = 12;
nTrajs_list =  baseNTrajPts * (startIdx : endIdx);


figure (55)

plot(nTrajs_list, errorResult.firstEuler(startIdx : endIdx),"-", Marker="o",LineWidth=1.5);
hold on
plot(nTrajs_list, errorResult.secondEuler(startIdx : endIdx), "-", Marker="<",LineWidth=1.5);
plot(nTrajs_list, errorResult.firstRk4(startIdx : endIdx), "-", Marker="s",LineWidth=1.5);
plot(nTrajs_list, errorResult.secondRk4(startIdx : endIdx), "-", Marker="v",LineWidth=1.5);
set(gca, 'YScale', 'log')
grid on
legend("1st-Euler", "2nd-Euler","1st-RK4","2nd-RK4",'Interpreter','latex','FontSize',10)
ylabel('Total Transcription Error $\eta_{total}$','Interpreter','latex','FontSize',12)
xlabel('$N$','Interpreter','latex','FontSize',12)


figure (66)

plot(nTrajs_list, timeResult.firstEuler(startIdx : endIdx), "-", Marker="o",LineWidth=1.5);
hold on
plot(nTrajs_list, timeResult.secondEuler(startIdx : endIdx), "-", Marker="<",LineWidth=1.5);
plot(nTrajs_list, timeResult.firstRk4(startIdx : endIdx), "-", Marker="s",LineWidth=1.5);
plot(nTrajs_list, timeResult.secondRk4(startIdx : endIdx), "-", Marker="v",LineWidth=1.5);
% set(gca, 'YScalte', 'log')
% ylim([0 0.8])
grid on
legend("1st-Euler", "2nd-Euler","1st-RK4","2nd-RK4",'Interpreter','latex','FontSize',10)
ylabel('NLP Solver Time $(s)$','Interpreter','latex','FontSize',12)
xlabel('$N$','Interpreter','latex','FontSize',12)


startIdx = 10;
nTrajs_list =  baseNTrajPts * (startIdx : endIdx);

figure (88)
grid on
plot(nTrajs_list, errorResult.firstEuler(startIdx :endIdx), Marker="o",LineWidth=1.5);
hold on
grid on
plot(nTrajs_list, errorResult.secondEuler(startIdx : endIdx), Marker="o",LineWidth=1.5);
plot(nTrajs_list, errorResult.firstRk4(startIdx : endIdx), Marker="o",LineWidth=1.5);
plot(nTrajs_list, errorResult.secondRk4(startIdx : endIdx), Marker="o",LineWidth=1.5);
ylim([0.3 1])
