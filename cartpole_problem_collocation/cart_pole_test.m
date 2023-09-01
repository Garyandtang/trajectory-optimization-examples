clc;clear;

%%% Setup casadi solver 
addpath(genpath("D:\software\casadi-windows-matlabR2016a-v3.5.5"))
import casadi.*
addpath("..\models\")
addpath("..\");
addpath("..\solvers\")

%%% setup dynamics config
config.dyn.g = 9.81;
config.dyn.l = 1;
config.dyn.m1 = 1;
config.dyn.m2 = 0.1;

%%% setup problem
problem = cart_pole_swing_up(config);

baseNTrajPts = 400;
% flag
config.flag.animationOn = true;
config.method.objAppro = "trapzoid_explict";

%%% solve problem
% % first euler
% config.method.dynamics = "first_order_euler";
% euler1Soln = directTranscriptionMethod(problem, config);
% 
% % second euler
% config.method.dynamics = "second_order_euler";
% euler2Soln = directTranscriptionMethod(problem, config);
% 
% % first rk4
% config.method.dynamics = "first_order_rk4";
% rk1Soln = directTranscriptionMethod(problem, config);

% second rk4
config.method.dynamics = "second_order_rk4";
rk2Soln = directTranscriptionMethod(problem, config);

%%% draw animation
soln = rk2Soln; % assign the soln to be drawnw 
P.plotFunc = @(t,z)( drawCartPole(t,z,config.dyn) );
P.speed = 0.7;
P.figNum = 102;
t = linspace(soln.tSoln(1),soln.tSoln(end),250);
z = soln.interp.x(t);
if config.flag.animationOn
   animate(t,z,P)
end

