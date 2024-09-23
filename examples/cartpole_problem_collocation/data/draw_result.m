%%%
% this script is to draw the result of 4 different shooting method for cart
% pole swing up problem.
%%%
clear; clc;
close all

load("cart_pole_result.mat")
line_width = 1.5;
startIdx = 2;
endIdx = 12;
nTrajs_list =  baseNTrajPts * (startIdx : endIdx);

figure (55)
plot(nTrajs_list, errorResult.cg1(startIdx : endIdx), Marker="o", LineWidth=line_width);
hold on
% set(gca, 'YScale', 'log')
plot(nTrajs_list, errorResult.cg2(startIdx : endIdx), Marker="^", LineWidth=line_width);
plot(nTrajs_list, errorResult.cg3(startIdx : endIdx), Marker="v", LineWidth=line_width);
plot(nTrajs_list, errorResult.ours(startIdx : endIdx), Marker="square", LineWidth=line_width);

legend("CG1", "CG2","CG3","Ours",'Interpreter','latex','FontSize',10)
ylabel('Total Transcription Error $\eta_{total}$','Interpreter','latex','FontSize',16)
xlabel('$N$','Interpreter','latex','FontSize',16)
grid on


figure (66)
plot(nTrajs_list, timeResult.cg1(startIdx : endIdx), Marker="o", LineWidth=line_width);
hold on
plot(nTrajs_list, timeResult.cg2(startIdx : endIdx), Marker="^", LineWidth=line_width);
plot(nTrajs_list, timeResult.cg3(startIdx : endIdx), Marker="v", LineWidth=line_width);
plot(nTrajs_list, timeResult.ours(startIdx : endIdx), Marker="square", LineWidth=line_width);
grid on
legend("CG1", "CG2","CG3","Ours",'Interpreter','latex','FontSize',10)
ylabel('NLP Solver Time $(s)$','Interpreter','latex','FontSize',16)
xlabel('$N$','Interpreter','latex','FontSize',16)

startIdx = 5;
endIdx = 12;
nTrajs_list =  baseNTrajPts * (startIdx : endIdx);

figure (88)
plot(nTrajs_list, errorResult.cg1(startIdx : endIdx), Marker="o", LineWidth=line_width);
hold on
plot(nTrajs_list, errorResult.cg2(startIdx : endIdx), Marker="^", LineWidth=line_width);
plot(nTrajs_list, errorResult.cg3(startIdx : endIdx), Marker="v", LineWidth=line_width);
plot(nTrajs_list, errorResult.ours(startIdx : endIdx), Marker="square", LineWidth=line_width);
xlim([baseNTrajPts*startIdx baseNTrajPts*endIdx])
grid on

% legend("CG1", "CG2","CG3","Ours",'Interpreter','latex','FontSize',10)
% ylabel('Total Transcription Error $\eta_{total}$','Interpreter','latex','FontSize',12)
% xlabel('$N$','Interpreter','latex','FontSize',12)