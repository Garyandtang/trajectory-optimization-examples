function traj = dirTrans_Trapz(config)
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
%   No of control input:
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
n = config.grid.nTrajPts;

% Create the initial guess:
guess.time = linspace(config.guess.time(1),config.guess.time(end),n);
guess.state = interp1(...
    config.guess.time', config.guess.state', guess.time')';
guess.control = interp1(...
    config.guess.time', config.guess.control', guess.time')';
% guess.state = zeros(size(guess.state));
% guess.control = zeros(size(guess.control));
% Create a list of all linear constraints, then add function handles:
[problem, pack] = buildConstraints(guess,config);
problem.objective = @(z)( costFunctionWrapper(z,pack) ); 
problem.nonlcon = @(z)( nonLinCon(z,pack, config) );

% Solve using fmincon:
[zSoln,fSoln,exitFlag] = fmincon(problem);

% Post-processing:
[t,x,u] = unPackDecVar(zSoln,pack);

traj.time = linspace(t(1),t(2),n);
traj.state = x;
traj.control = u;
traj.objVal = fSoln;
traj.exitFlag = exitFlag;
traj.interp.state = @(tt)( interp1(traj.time',traj.state',tt')' );
traj.interp.control = @(tt)( interp1(traj.time',traj.control',tt')' );

end

function [C, Ceq] = nonLinCon(z,pack,config)
%
% This function enforces the dynamics of the cart-pole system
%

% Based on the Trapazoid Method for discretization, as defined in Bett's book,
% chapter 4.

n = pack.nState(2);
[t,x,u] = unPackDecVar(z,pack);
dt = (t(2)-t(1))/(n-1);

% Evaluate the dynamics at each collocation point
dx = cartPoleDynamics(x,[u],config.dyn);  

% Trapazoid rule:
idxLow = 1:(n-1);
idxUpp = 2:n;
intStateTrap = 0.5*dt*(dx(:,idxLow) + dx(:,idxUpp));
intStateCol = x(:,idxUpp)-x(:,idxLow);

% Defect constraint:
defect = intStateTrap - intStateCol;
% 
% % user-defined boundary constraints:
% [bndIneq, bndEq] = boundaryConstraint(t,x(:,1),x(:,end),config.userData);

C = [];%bndIneq;
Ceq = [reshape(defect,numel(defect),1)];% bndEq];


end


function cost = costFunctionWrapper(z,pack)

[t,x,u] = unPackDecVar(z,pack);

% Trapazoid rule to integrate cost function:
tt = linspace(t(1),t(end),pack.nState(2));  %Time vector
dc = costFunction(tt,x,u);
nTime = size(x,2);
dt = (t(2)-t(1))/(nTime-1);
w = ones(nTime,1); 
w([1,end]) = 0.5;
cost = dt*dc*w;

end

