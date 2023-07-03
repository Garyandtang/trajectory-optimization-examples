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
%                               q(1) = 1,   dq(1)= 0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc;clear;

%%% Setup casadi solver 
addpath(genpath("D:\software\casadi-windows-matlabR2016a-v3.5.5"))
import casadi.*
addpath("..\models\")
addpath("..\")
addpath('utils');

% setup dynamics config 
config.m1 = 1;
config.m2 = 1;
config.g = 9.81;
config.l1 = 0.5;
config.l2 = 0.5;

%%% setup problem
problem = acrobot_swing_up(config);

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
    % Interpolate the solution on a uniform grid for plotting and animation:
    tGrid = soln.tSoln;
    t = linspace(tGrid(1),tGrid(end),1000);
    q = soln.interp.q(t);
    dq = soln.interp.dq(t);
    z = [q;dq];
    u = soln.interp.u(t);
    % Animate the results:
    A.plotFunc = @(t,z)( drawAcrobot(t,z,config) );
    A.speed = 0.25;
    A.figNum = 101;
    animate(t,z,A)
    
    % Plot the results:
    figure(1337); clf; plotAcrobot(t,z,u,config);
end






