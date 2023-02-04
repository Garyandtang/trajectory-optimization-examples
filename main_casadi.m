% main_casadi.m
%
% This script runs general direct collocation methods for trajectory
% optimization (trapazoid) with casadi framework
% it uses augumented state to cast second-order system into a first-order 
% form.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem description
%
% Cartpole swing up: 
%   State:
%       z = [4, 1] = state vector = [x;q;dx;dq]
%   
%   Control:
%       u = [1, 1] = actuation vector = F = force on cart
%
%   objective:
%       L = u^2

%
%   Parameters:
%       .m1 = cart mass
%       .m2 = pendulum point-mass
%       .g = gravity
%       .l = length of the pendulum
%     
%   System dynamics:
%       dz = f(z, u)
%       dz(1) = z(3)  % first order
%       dz(2) = z(4)  % first order
%       dz(3) = ddx
%       dz(4) = ddq
%       ddx = l*m2*sin(q)*dq2^2 + u + m2*g*cos(q)*sin(q)/(m1+m2*(1 - cos(q)^2))
%       ddq = -(l*m2*cos(q)*sin(q)*dq2^2 + u*cos(q)+(m1+m2)*g*sin(q))/(l*m1+l*m2*(1 - cos(q)^2))
%
%
%   Constraints:
%       -49 <= u <= 49      bounded control input 
%       -1 <= x <= 1        bounded position
%       0 <= tf - t0 <= x   bounded duration 
%       z_init = [0,0,0,0]  init state constraint
%       z_f = [0,pi,0,0]    final state constraint
%
%   No of grid points:
%       N = 15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direction collocation formulation:
% Trapezoid method:
%   Assume dynamics and control are linear between grid points 
%   
%   Trapezoidal rule (left and right Rieman sums):
%       \int_a^b f(x)dx = (b-a)*0.5*(f(a) + f(b))
%
%   Decision variable:
%       w = [2 + 4*15 + 1*15, 1] = [t0; tf; X; U];
%   
%   Bounded contraints (as above):
%       
%   defect constraints (by trapezoidal rule):
%       xk+1 = xk + hk/2(fk+1 + fk)
%
%   objective:
%       L = u^2
%       J = \int u(t)^2dt   continuous form
%       J = \sum hk/2(uk + uk+1) discreted form by trapezoid quadrature
%
%   Parameters:
%       .m1 = cart mass
%       .m2 = pendulum point-mass
%       .g = gravity
%       .l = length of the pendulum
%     
%   System dynamics:
%       dx0 = (1 - x2^2)*x1 - x2 + u
%       dx1 = x0
%
%   Constraints:
%       bounded control input 
%       bounded position
%       bounded time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear;
%%% Setup casadi solver 
addpath(genpath("D:\software\casadi-windows-matlabR2016a-v3.5.5"))
import casadi.*

% flag
singleTimeDV = 1;

% Number of points for initialization:
config.grid.nTrajPts = 20;

% Physical parameters for dynamics
m1 = 2.0; config.dyn.m1 = m1;   %cart mass
m2 = 0.5; config.dyn.m2 = m2;   %pendulum mass
config.dyn.g = 9.81;
config.dyn.l = 1;
g = config.dyn.g;
l = config.dyn.l;

% Compute an initial guess at a trajectory:
config.guess = computeGuess(config);

% Bounds:
config.bounds = computeBounds(config);
bounds = config.bounds;

%%% Declare model variables
x = MX.sym('x');    % cart position
q = MX.sym('q');    % pole angle
dx = MX.sym('dx');  % cart velocity
dq = MX.sym('dq');  % pole angular rate
z = [x; q; dx; dq]; % state vector of cartpole
u = MX.sym('u');    % control input: force on cart

%%% Model equations: continuous system dynamics
ddx = (l*m2*sin(q)*dq^2 + u + m2*g*cos(q)*sin(q))/(m1+m2*(1 - cos(q)^2));
ddq = -(l*m2*cos(q)*sin(q)*dq^2 + u*cos(q)+(m1+m2)*g*sin(q))...
        /(l*m1+l*m2*(1 - cos(q)^2));
dz = [z(3);
      z(4);
      ddx;
      ddq];

f = Function('f', {z, u}, {dz}); % CT dynamics

%%% Formulate the NLP
% decistion variable: w = [1*15 + 4*15 + 1*15, 1] = [T; X; U];
numDecVar = (1 + size(z, 1) + size(u, 1)) * config.grid.nTrajPts;
if singleTimeDV
    w0 = [config.guess.time(1);
          config.guess.time(end);
          reshape(config.guess.state,numel(config.guess.state), 1);
          reshape(config.guess.control,numel(config.guess.control), 1)];
else
    w0 = [reshape(config.guess.time, numel(config.guess.time), 1);
          reshape(config.guess.state,numel(config.guess.state), 1);
          reshape(config.guess.control,numel(config.guess.control), 1)];
end

% declare decision variables and constraints
decVar_T = [];  % all declare decision on time
decVar_Z = [];  % all declare decision on state
decVar_U = [];  % all declare decision on control
g = [];         % all nonlinear constraints function for IPOPT
lbt = [];       % all lower bounds on time
ubt = [];       % all upper bounds on time
lbz = [];       % all lower bounds on state
ubz = [];       % all upper bounds on state
lbu = [];       % all lower bounds on control
ubu = [];       % all upper bounds on control
lbg = [];       % all lower bounds on g
ubg = [];       % all upper bounds on g

% decision variable on T and bounded time constaints
if singleTimeDV
    decVar_T = MX.sym('T', 2);
    lbt = [bounds.initialTime.lower; bounds.finalTime.lower];
    ubt = [bounds.initialTime.upper; bounds.finalTime.upper];
else
    decVar_T = [decVar_T; MX.sym('T', config.grid.nTrajPts)];
    lbt = [lbt; bounds.initialTime.lower];
    ubt = [ubt; bounds.initialTime.upper];
    for i = 1 : (config.grid.nTrajPts - 2)
        lbt = [lbt; bounds.initialTime.lower];
        ubt = [ubt; bounds.finalTime.upper];
    end
    lbt = [lbt; bounds.finalTime.lower];
    ubt = [ubt; bounds.finalTime.upper];
end

% decision variable initial z and u and bounded contraints
zk = MX.sym(['X_' num2str(1)], 4);
uk = MX.sym(['U_' num2str(1)]);
dzk = f(zk, uk);
decVar_Z = [decVar_Z; zk];     
decVar_U = [decVar_U; uk];      
lbu = [lbu; bounds.control.lower];   % concatenated control lower bound
ubu = [ubu; bounds.control.upper];   % concatenated control upper bound
lbz = [lbz; bounds.state.lower];     % concatenated state lower bound 
ubz = [ubz; bounds.state.upper];     % concatenated state upper bound

g = [g;zk];                     % add init state bounds to g
lbg = [lbg; zeros(4,1)];        % add lower bound for init state
ubg = [ubg; zeros(4,1)];        % add upper bound for init state



J = 0;
% decision variable on path z and u and bounded contraints + defect
% constraints
dt = (decVar_T(2) - decVar_T(1)) / (config.grid.nTrajPts - 1);
for k = 2 : config.grid.nTrajPts
    dzk_prev = dzk;
    zk_prev = zk;
    uk_prev = uk;
    dzk = f(zk, uk);
    zk = MX.sym(['X_' num2str(k)], 4);
    uk = MX.sym(['U_' num2str(k)]);
    decVar_U = [decVar_U; uk];
    decVar_Z = [decVar_Z; zk];
    lbu = [lbu; bounds.control.lower];
    ubu = [ubu; bounds.control.upper];
    lbz = [lbz; bounds.state.lower];
    ubz = [ubz; bounds.state.upper];
    % calculate cost function by tapezoidal rule
    J = J + dt*(uk_prev^2 + uk^2)/2;

    % calcuate defect constraints besed on trapezoidal rule
    g1 = zk - zk_prev - dt/2*(dzk + dzk_prev);
    g = [g; g1];
    lbg = [lbg; zeros(4, 1)];
    ubg = [ubg; zeros(4,1)];

end
% add final state bounded
g = [g; zk];
lbg = [lbg;-0.7;pi;0;0];
ubg = [ubg;-0.7;pi;0;0];

% concatenate all decVar
w = [decVar_T; decVar_Z; decVar_U];
lbw = [lbt; lbz; lbu];
ubw = [ubt; ubz; ubu];

% Create an NLP solver
prob = struct('f', J, 'x', w, 'g', g);
solver = nlpsol('solver', 'ipopt', prob);

% Solve the NLP
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);

dim.nTime = [1, size(decVar_T, 1)];
dim.nState = [4, config.grid.nTrajPts];
dim.nControl = [1, config.grid.nTrajPts];
% Post-processing:
[t,x,u] = unPackDecVar(w_opt,dim);
traj.time = linspace(t(1),t(2),config.grid.nTrajPts);
traj.state = x;
traj.control = u;
traj.objVal = 1;
traj.exitFlag = 1;
traj.interp.state = @(tt)( interp1(traj.time',traj.state',tt')' );
traj.interp.control = @(tt)( interp1(traj.time',traj.control',tt')' );

P.plotFunc = @(t,z)( drawCartPole(t,z,config.dyn) );
P.speed = 0.7;
P.figNum = 102;
t = linspace(traj.time(1),traj.time(end),250);
z = traj.interp.state(t);
animate(t,z,P)

% Plot the results:
figure(101); clf; plotTraj(traj,config);