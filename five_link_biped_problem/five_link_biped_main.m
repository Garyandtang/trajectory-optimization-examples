% block_move_main.m
%
% This script is used to demo the result of block move example
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Continous time trajectory formuation
%
% model:
%   see: blockMoveModel.m
%
% problem:
%       min_{u(t),q(t),dq(t)}  \int_0^1 u(t)^2 dt
%           s.t.                ddq(t) = u(t)
%                               q(0) = 0,   dq(1)=0,\\
%                               q(1) = 1,trajhandle   dq(1)= 0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc;clear;

%%% Setup casadi solver 
addpath(genpath("D:\software\casadi-windows-matlabR2016a-v3.5.5"))
import casadi.*
addpath("..\models\")
addpath("..\")
addpath('utils');

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
config.method.dynamics = "first_order_euler";
euler1Soln = directTranscriptionMethod(problem, config);
% 
% % second euler
% config.method.dynamics = "second_order_euler";
% euler2Soln = directTranscriptionMethod(problem,config);
% 
% % first rk4
% config.method.dynamics = "first_order_rk4";
% rk1Soln = directTranscriptionMethod(problem,  config);

% second rk4
config.method.dynamics = "second_order_rk4";
rk2Soln = directTranscriptionMethod(problem,  config);

if config.flag.animationOn
    soln = rk2Soln;
    q = soln.qSoln;
    t = soln.tSoln;
    % Interpolate the solution on a uniform grid for plotting and animation:
    Anim.speed = 0.25;
    Anim.plotFunc = @(t,q)( drawRobot(q,config) );
    Anim.verbose = true;
    animate(t,q,Anim);

end






