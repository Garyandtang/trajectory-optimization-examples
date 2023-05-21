function model = cartPoleModel(config)
import casadi.*

if (~exist('config'))
    warning("Config does not include carport dynamic, use default one.")
    g = 9.81;
    l = 1;
    m1 = 1;
    m2 = 1;
else
    % get model configuration
    g = config.dyn.g;
    l = config.dyn.l;
    m1 = config.dyn.m1;
    m2 = config.dyn.m2;
end
%%% Declare model variables
q1 = MX.sym('q1');    % cart position
q2 = MX.sym('q2');    % pole angle
dq1 = MX.sym('dq1');  % cart velocity
dq2 = MX.sym('dq2');  % pole angular rate
q = [q1; q2];         % configuration vector
dq = [dq1; dq2];      % first time-deriviative of configuration        
u = MX.sym('u');      % control input: force on car
dt = MX.sym('dt');    % 
z = [q;dq];
%%% Model equations: continuous system dynamics
ddq1 = (l*m2*sin(q2)*dq2^2 + u + m2*g*cos(q2)*sin(q2))/(m1+m2*(1 - cos(q2)^2));
ddq2 = -(l*m2*cos(q2)*sin(q2)*dq2^2 + u*cos(q2)+(m1+m2)*g*sin(q2))...
        /(l*m1+l*m2*(1 - cos(q2)^2));
ddq = [ddq1; ddq2];
dz = [dq; ddq];

f = Function('f', {z, u}, {dz}); % CT dynamics
f2 = Function('f2', {q,dq, u}, {ddq});

% CT dynamic jacobian matrix
f_x = jacobian(dz, z);   % 4 * 4
f_u = jacobian(dz, u);   % 4 * 1

model.CTDynamic.evaluation.firstOrder = f;
model.CTDynamic.evaluation.secondOrder = f2;
% model.CTDynamic.jacobian.f_x = f_x;
% model.CTDynamic.jacobian.f_u = f_u;

%%% Model equations: discrete system dynamics
%%% Euler Method
% Evaluation
zNext = z + dt*dz;
f_Euler = Function('f', {z, u, dt}, {zNext}); % DT dynamics

% jacobian matrix
f_x_dt = eye(size(f_x)) + dt * f_x;
f_u_dt = dt * f_u;
A_hat = Function('A_hat', {z, u, dt}, {f_x_dt});
B_hat = Function('B_linear', {z, u, dt}, {f_u_dt});

model.DTDynamic.Euler.evaluation = f_Euler;
model.DTDynamic.Euler.jacobain.f_x = A_hat;
model.DTDynamic.Euler.jacobain.f_u = B_hat;

model.dim.nState = 4;
model.dim.nControl = 1;
model.dim.nConfig = 2;

end