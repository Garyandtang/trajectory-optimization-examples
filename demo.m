clc;clear;

%%% Setup casadi solver 
import casadi.*
addpath('examples\five_link_biped_problem\utils\');
addpath('examples\five_link_biped_problem\');
% setup dynamics config 
config.m1 = 3.2;
config.m2 = 6.8;
config.m3 = 20;
config.m4 = 6.8;
config.m5 = 3.2;
config.I1 = 0.93;
config.I2 = 1.0800;
config.I3 = 2.2200;
config.I4 = 1.0800;
config.I5 = 0.9300;
config.l1 = 0.4000;
config.l2 = 0.4000;
config.l3 = 0.6250;
config.l4 = 0.4000;
config.l5 = 0.4000;
config.c1 = 0.1280;
config.c2 = 0.1630;
config.c3 = 0.2000;
config.c4 = 0.1630;
config.c5 = 0.1280;
config.g = 9.8100;
config.stepLength = 0.5000;
config.stepTime = 0.7000;

%%% setup problem
problem = five_link_biped_walk(config);

% flag
config.flag.animationOn = true;
config.method.objAppro = "trapzoid_explict";


% first euler
problem.grid.nTrajPts = 100;
config.method.dynamics = "second_order_trapzoidal";
soln = directTranscriptionMethod(problem, config);

if config.flag.animationOn
    q = soln.qSoln;
    t = soln.tSoln;
    % Interpolate the solution on a uniform grid for plotting and animation:
    Anim.speed = 0.25;
    Anim.plotFunc = @(t,q)( drawRobot(q,config) );
    Anim.verbose = true;
    animate(t,q,Anim);

end
