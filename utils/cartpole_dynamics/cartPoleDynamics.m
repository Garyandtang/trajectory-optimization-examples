function dz = cartPoleDynamics(z,u,p)
% dz = cartPoleDynamics(z,u,p)
%
% This function computes the dynamics for a simple cart-pole: a point-mass
% cart that travels on a friction-less horizontal track. A pendulum hangs
% from the cart. There are two actuators: a horizontal force that pushes on
% the cart, and a torque acting on the base of the pendulum.
%
% INPUTS:
%   z = [4, n] = state vector = [x;q;dx;dq]
%   u = [1, n] = actuation vector = F = force on cart
%   p = struct of parameters:
%       .m1 = cart mass
%       .m2 = pendulum point-mass
%       .g = gravity
%       .l = length of the pendulum
%
% OUTPUTS:
%   dz = [4, n] = derivative of the state vector
%

% x = z(1,:);
q = z(2,:);
dx = z(3,:);
dq = z(4,:);

F = u(1,:);
T = zeros(1,length(F));

% [ddx,ddq] = autoGen_cartPoleDynamics(q,dq,F,T,p.m1,p.m2,p.g,p.l);
l = p.l;
g = p.g;
m2 = p.m2;
m1 = p.m1;
ddx = zeros(1, size(z, 2));
ddq = zeros(1, size(z, 2));
for i = 1 : size(z,2)
    q = z(2, i);
    dq2 = z(4, i);
    u = F(1, i);
    ddx(1,i) = l*m2*sin(q)*dq2^2 + u + m2*g*cos(q)*sin(q)/(m1+m2*(1 - cos(q)^2));
    ddq(1,i) = -(l*m2*cos(q)*sin(q)*dq2^2 + u*cos(q)+(m1+m2)*g*sin(q))/(l*m1+l*m2*(1 - cos(q)^2));
end
dz = [dx;dq;ddx;ddq];

end


