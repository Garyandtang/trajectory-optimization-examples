function ddq_ = cartPoleDynamicsCasadi(config,q_, dq_, u_)
addpath(genpath("D:\software\casadi-windows-matlabR2016a-v3.5.5"))
import casadi.*
g = config.dyn.g;
l = config.dyn.l;
m1 = config.dyn.m1;
m2 = config.dyn.m2;

%%% Declare model variables
q1 = MX.sym('q1');    % cart position
q2 = MX.sym('q2');    % pole angle
dq1 = MX.sym('dq1');  % cart velocity
dq2 = MX.sym('dq2');  % pole angular rate
q = [q1; q2];         % configuration vector
dq = [dq1; dq2];      % first time-deriviative of configuration        
u = MX.sym('u');    % control input: force on cart

%%% Model equations: continuous system dynamics
ddq1 = (l*m2*sin(q2)*dq2^2 + u + m2*g*cos(q2)*sin(q2))/(m1+m2*(1 - cos(q2)^2));
ddq2 = -(l*m2*cos(q2)*sin(q2)*dq2^2 + u*cos(q2)+(m1+m2)*g*sin(q2))...
        /(l*m1+l*m2*(1 - cos(q2)^2));
ddq = [ddq1; ddq2];
f = Function('f', {q,dq, u}, {ddq});

%%% get the results
ddq_ = f(q_,dq_,u_);
ddq_ = full(ddq_);
end