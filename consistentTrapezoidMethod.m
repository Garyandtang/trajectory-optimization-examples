function soln = consistentTrapezoidMethod(config)
% consistentTrapezoidMethod.m
%
% This script runs consistent direct collocation methods for trajectory
% optimization (trapazoid) with casadi framework
%
% implement as functio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem description
%
% Cartpole swing up: 
%   configuration:
%       q = [2, 1] = configuration vector = [q1,q2]
%       q1: postion of cart (m)
%       q2: angle of pole (rad)
%
%   first-order derivative of configuration
%       dq = [2, 1] = [dq1,dq2]
%       dq1: velocity of cart (m/s)
%       dq2: angular velocity of pole (rad/s)
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
%       ddq = g(q, dq, u)
%       ddq1 = (l*m2*sin(q2)*dq2^2 + u + m2*g*cos(q2)*sin(q2))/(m1+m2*(1 - cos(q2)^2));
%       ddq2 = -(l*m2*cos(q2)*sin(q2)*dq2^2 + u*cos(q2)+(m1+m2)*g*sin(q2)) /(l*m1+l*m2*(1 - cos(q2)^2));
%
%   Constraints:
%       -umax <= u <= umax      bounded control input 
%       -q1max <= x <= q1max    bounded position
%       0 <= tf - t0 <= x       bounded duration 
%       z_init = [0,0,0,0]      init state constraint
%       z_f = [0,pi,0,0]        final state constraint
%
%   No of grid points:
%       N = 15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consisent Direction collocation formulation:
% Trapezoid method:
%   See paper Section IV-C
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setup casadi solver 
addpath(genpath("D:\software\casadi-windows-matlabR2016a-v3.5.5"))
import casadi.*

% flag
singleTimeDV = 1;
plotError = 0;

% get cartpole model
model = cartPoleModel(config);
nState = model.dim.nState;
nControl = model.dim.nControl;
nConfig = model.dim.nConfig;
f = model.CTDynamic.evaluation.secondOrder;

% get problem configuration
bounds = config.bounds;
initState = config.initState;
finalState = config.finalState;

%%% Formulate the NLP
% decistion variable: w = [1*15 + 4*15 + 1*15, 1] = [T; X; U];
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
qk = MX.sym(['Q_' num2str(1)], 2);
vk = MX.sym(['V_' num2str(1)], 2);
uk = MX.sym(['U_' num2str(1)]);
ddqk = f(qk, vk, uk);
decVar_Z = [decVar_Z; qk; vk];     
decVar_U = [decVar_U; uk];      
lbu = [lbu; bounds.control.lower];   % concatenated control lower bound
ubu = [ubu; bounds.control.upper];   % concatenated control upper bound
lbz = [lbz; bounds.state.lower];     % concatenated state lower bound 
ubz = [ubz; bounds.state.upper];     % concatenated state upper bound

g = [g;qk;vk];                     % add init state bounds to g
lbg = [lbg; initState];        % add lower bound for init state
ubg = [ubg; initState];        % add upper bound for init state



J = 0;
% decision variable on path z and u and bounded contraints + defect
% constraints
dt = (decVar_T(2) - decVar_T(1)) / (config.grid.nTrajPts - 1);
for k = 2 : config.grid.nTrajPts
    qk_prev = qk;
    vk_prev = vk;
    ddqk_prev = ddqk;
    uk_prev = uk;
    qk = MX.sym(['Q_' num2str(1)], 2);
    vk = MX.sym(['V_' num2str(1)], 2);
    uk = MX.sym(['U_' num2str(k)]);
    ddqk = f(qk, vk, uk);
    decVar_Z = [decVar_Z; qk; vk];
    decVar_U = [decVar_U; uk];
    lbu = [lbu; bounds.control.lower];
    ubu = [ubu; bounds.control.upper];
    lbz = [lbz; bounds.state.lower];
    ubz = [ubz; bounds.state.upper];
    % calculate cost function based on approximated u(t)
    J = J + dt*(uk_prev^2 + uk^2 + uk_prev * uk)/3;
%     J = J + dt*(uk_prev^2 + uk^2 )/2;

    % calcuate defect constraints besed on consistent trapezoidal rule
    g1 = vk - vk_prev -dt/2*(ddqk + ddqk_prev);
    g2 = qk - qk_prev -vk_prev*dt - dt^2*(2*ddqk_prev + ddqk)/6;
    g = [g; g1;g2];
    lbg = [lbg; zeros(4, 1)];
    ubg = [ubg; zeros(4,1)];

end
% add final state bounded
g = [g; qk;vk];
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
soln.objVal = full(sol.f);
soln.exitFlag = 1;
% traj.interp.state = @(tt)( interp1(traj.time',traj.state',tt')' );
soln.interp.control = @(tt)( interp1(soln.time',soln.control',tt')' );

% use piecewise quadratic interpolation first-derivative of configuration
% use piecewise cubic interpolation for configuration
soln.ddConfig = full(f(soln.configuration, soln.dConfiguration, soln.control));
soln.interp.dConfig = @(tt)( bSpline2(soln.time, soln.dConfiguration, ...
    soln.ddConfig, tt));
soln.interp.config = @(tt) (bSpline3(soln.time, soln.configuration,...
    soln.dConfiguration, soln.ddConfig, tt));
soln.interp.state = @(tt) ([soln.interp.config(tt); soln.interp.dConfig(tt)]);
soln.interp.ddConfig = @(tt) (interp1(soln.time',soln.ddConfig',tt)');

% interpolation for checking collocation constraint along the trjeatory
%  collocation constraint = dq(t) - v(t), here v = dq, so it is zeros
soln.interp.collCst = @(t)(0.*t.*zeros(2,1));

% the result here is zero is that we asume u(t) is a piecewise linear
% function and we derive the approximated objective function based on
% u(t), and the objective function just depends on u(t), hence the objCst
% is zero. In general objective function \int w(u,x,t) dt, for example
% w = xQx + uRu. Then objCst will not be zero
soln.interp.objCst = @(t)(t*0);

% use romberg quadrature to estimate the absolute dynamic error
absColErr = @(t)(abs(soln.interp.collCst(t)));
nSegment = config.grid.nTrajPts-1;
nState = size(soln.configuration,1);
quadTol = 1e-12;   %Compute quadrature to this tolerance  
soln.info.dynError = zeros(nState,nSegment);
for i=1:nSegment
    soln.info.dynError(:,i) = rombergQuadrature(absColErr,soln.time([i,i+1]),quadTol);
end
soln.info.maxError = max(max(soln.info.dynError));

% use romberg quadrature to estimate the absolute objective error
absObjErr = @(t)(abs(soln.interp.objCst(t)));
nSegment = config.grid.nTrajPts-1;
nObj = 1;
quadTol = 1e-12;   %Compute quadrature to this tolerance  
soln.info.objError = zeros(nObj,nSegment);
for i=1:nSegment
    soln.info.objError(:,i) = rombergQuadrature(absObjErr,soln.time([i,i+1]),quadTol);
end
soln.info.maxObjError = max(max(soln.info.objError));


P.plotFunc = @(t,z)( drawCartPole(t,z,config.dyn) );
P.speed = 0.7;
P.figNum = 102;
t = linspace(soln.time(1),soln.time(end),250);
z = soln.interp.state(t);
if config.flag.animationOn
   animate(t,z,P)
end


% % Plot the results:
% figure(101); clf; plotTraj(soln,config);

%%%% Show the error in the collocation constraint between grid points:
%
if 0
    % Then we can plot an estimate of the error along the trajectory
    figure(5); clf;
    
    % NOTE: the following commands have only been implemented for the direct
    % collocation(trapezoid, hermiteSimpson) methods, and will not work for
    % chebyshev or rungeKutta methods.
    cc = soln.interp.collCst(t);
    
    subplot(2,2,1);
    plot(t,cc(1,:))
    title('Collocation Error:   dx/dt - f(t,x,u)')
    ylabel('d/dt cart position')
    
    subplot(2,2,3);
    plot(t,cc(2,:))
    xlabel('time')
    ylabel('d/dt pole angle')
    
    idx = 1:length(soln.info.dynError);
    subplot(2,2,2); hold on;
    plot(idx,soln.info.dynError(1,:),'ko');
    title('State Error')
    ylabel('cart position')
    
    subplot(2,2,4); hold on;
    plot(idx,soln.info.dynError(2,:),'ko');
    xlabel('segment index')
    ylabel('pole angle');
end
end