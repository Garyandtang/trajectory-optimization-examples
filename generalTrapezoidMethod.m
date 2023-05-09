function soln = generalTrapezoidMethod(config)
% main_casadi.m
%
% This script runs general direct collocation methods for trajectory
% optimization (trapazoid) with casadi framework
% it uses augumented state to cast second-order system into a first-order 
% form.
%
% implement as functio
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
%   See paper Section-III
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setup casadi solver 
addpath(genpath("D:\software\casadi-windows-matlabR2016a-v3.5.5"))
import casadi.*
addpath("models\")
% flag
singleTimeDV = 1;
plotError = 0;

% get cartpole model
model = cartPoleModel(config);
nState = model.dim.nState;
nControl = model.dim.nControl;
nConfig = model.dim.nConfig;
f = model.CTDynamic.evaluation.firstOrder;
f2 = model.CTDynamic.evaluation.secondOrder;

% get problem configuration
bounds = config.bounds;
initState = config.initState;
finalState = config.finalState;

%%% Formulate the NLP
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

g = [g; zk];                     % add init state bounds to g
lbg = [lbg; initState];        % add lower bound for init state
ubg = [ubg; initState];        % add upper bound for init state



J = 0;
% decision variable on path z and u and bounded contraints + defect
% constraints
dt = (decVar_T(2) - decVar_T(1)) / (config.grid.nTrajPts - 1);
for k = 2 : config.grid.nTrajPts
    zk_prev = zk;
    uk_prev = uk;
    dzk_prev = dzk;
    zk = MX.sym(['X_' num2str(k)], 4);
    uk = MX.sym(['U_' num2str(k)]);
    dzk = f(zk, uk);
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
lbg = [lbg;finalState];
ubg = [ubg;finalState];

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
[t,soln.state,soln.control] = unPackDecVar(w_opt,dim);
soln.time = linspace(t(1),t(2),config.grid.nTrajPts);
soln.configuration = soln.state(1:2,:);
soln.dConfiguration = soln.state(3:4,:);
tSoln = soln.time;
xSoln = soln.state;
uSoln = soln.control;
soln.ddConfiguration = full(f2(soln.configuration, soln.dConfiguration, soln.control));
soln.dJ = soln.control.*soln.control;
soln.objVal = full(sol.f);

soln.exitFlag = 1;
soln.interp.state = @(tt)( interp1(tSoln',xSoln',tt')' );
soln.interp.control = @(tt)( interp1(tSoln',uSoln',tt')' );

% use piecewise quadratic interpolation for the state
fSoln = full(f(xSoln, uSoln));
soln.interp.configuration = @(t)( bSpline2(tSoln,xSoln(1:2,:),fSoln(1:2,:),t) );
soln.interp.dConfiguration = @(t)( bSpline2(tSoln,xSoln(3:4,:),fSoln(3:4,:),t) );
soln.interp.ddConfiguration = @(t)(interp1(tSoln',fSoln(3:4,:)',t)');
soln.interp.state = @(t)( [soln.interp.configuration(t);soln.interp.dConfiguration(t)] );


soln.interp.collCst = @(t)(...
     soln.interp.dConfiguration(t) - interp1(tSoln',fSoln(1:2,:)',t)');

soln.interp.objCst = @(t)(interp1(tSoln',(uSoln.*uSoln)',t) ...
    - soln.interp.control(t).*soln.interp.control(t));

% use romberg quadrature to estimate the absolute dynamic error
absColErr = @(t)(abs(soln.interp.collCst(t)));
nSegment = config.grid.nTrajPts-1;
nState = size(nConfig,1);
quadTol = 1e-12;   %Compute quadrature to this tolerance  
soln.info.dynError = zeros(nConfig,nSegment);
for i=1:nSegment
    soln.info.dynError(:,i) = rombergQuadrature(absColErr,tSoln([i,i+1]),quadTol);
end
soln.info.maxdynError = max(max(soln.info.dynError));


% use romberg quadrature to estimate the absolute objective error
absObjErr = @(t)(abs(soln.interp.objCst(t)));
nSegment = config.grid.nTrajPts-1;
nObj = 1;
quadTol = 1e-12;   %Compute quadrature to this tolerance  
soln.info.objError = zeros(nObj,nSegment);
for i=1:nSegment
    soln.info.objError(:,i) = rombergQuadrature(absObjErr,tSoln([i,i+1]),quadTol);
end
soln.info.maxObjError = max(max(soln.info.objError));


% plot the animation
P.plotFunc = @(t,z)( drawCartPole(t,z,config.dyn) );
P.speed = 0.7;
P.figNum = 102;
t = linspace(tSoln(1),tSoln(end),100000);
z = soln.interp.state(t);
if config.flag.animationOn
   animate(t,z,P)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end