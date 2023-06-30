function model = quadrotor1d(config)
import casadi.*

% get config prameters
m = config.m;
g = config.g;

%%% Declare model variables
q = MX.sym('q');    % block postion
dq = MX.sym('dq');  % block velocity
u = MX.sym('u');    % control input: force on block (unit mass)
z = [q;dq];         % first-order state vector

%%% Model equation: Continuous time system dynamics
ddq = u/m - g;
dz = [dq;ddq];

f = Function('f', {z, u}, {dz});
f2 = Function('f2', {q, dq, u}, {ddq});


model.CTDynamic.evaluation.firstOrder = f;
model.CTDynamic.evaluation.secondOrder = f2;

%%% todo: jacobian and discrete time system dynamic


model.dim.nState = size(z, 1);
model.dim.nControl = size(u, 1);
model.dim.nConfig = size(q, 1);


end