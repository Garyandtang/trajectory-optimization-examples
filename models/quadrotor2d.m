function model = quadrotor2d(config)
import casadi.*
if (~exist("config"))
    m = 0.18;               % Gravitational acceleration (m/s^2)
    g = 9.81;               % Mass (kg)
    Ixx = 0.0025;           % Mass moment of inertia (kg*m^2)
    l = 0.086;              % Arm length (m)
    d = l / sqrt(2);
else
    % get config prameters
    m = config.m;
    g = config.g;
    d = config.l / sqrt(2);
    Ixx = config.Ixx;
end

%%% Declare model variables
q = MX.sym('q');    % block postion
dq = MX.sym('dq');  % block velocity
u = MX.sym('u');    % control input: force on block (unit mass)
z = [q;dq];         % first-order state vector

q1 = MX.sym('q1');      % x position
dq1 = MX.sym('dq1');    % x velocity
q2 = MX.sym('q2');      % z position
dq2 = MX.sym('dq2');    % z velocity
q3 = MX.sym('q3');      % theta
dq3 = MX.sym('dq3');    % angle rate
u1 = MX.sym('u1');      % thrust of left motor
u2 = MX.sym('u2');      % thrust of right motor

q = [q1;q2;q3];
dq = [dq1;dq2;dq3];
u =  [u1;u2];
x = [q;dq];


%%% Model equation: Continuous time system dynamics
ddq1 = sin(q3) * (u1 + u2) / m;
ddq2 = cos(q3) * (u1 + u2) / m - g;
ddq3 = (u2 - u1)*d / Ixx;
ddq = [ddq1; ddq2; ddq3];
dx = [dq;ddq];

f = Function('f', {x, u}, {dx});
f2 = Function('f2', {q, dq, u}, {ddq});


model.CTDynamic.evaluation.firstOrder = f;
model.CTDynamic.evaluation.secondOrder = f2;

%%% todo: jacobian and discrete time system dynamic


model.dim.nState = size(x, 1);
model.dim.nControl = size(u, 1);
model.dim.nConfig = size(q, 1);


end