function problem = acrobot_swing_up(config)
addpath("..\models\")
%%% setup dynamics
model = acrobot(config);
problem.model = model;
nState = model.dim.nState;
nControl = model.dim.nControl;

%%% setup cost
problem.cost.type = "LQR";
problem.cost.R = 1 * eye(nControl);
problem.cost.Q = 0 * eye(nState);

%%% setup constraints
% boundaries constraints
problem.constraints.boundarys.initState = [0;0;0;0];
problem.constraints.boundarys.finalState = [pi;pi;0;0];
problem.constraints.boundarys.finalTime = 2;

% bound constraints
problem.constraints.bounds.control.lower = -20*ones(nControl,1);
problem.constraints.bounds.control.upper = 20*ones(nControl,1);
problem.constraints.bounds.state.lower = [-2*pi; -2*pi; -inf(2,1)];
problem.constraints.bounds.state.upper = [ 2*pi;  2*pi;  inf(2,1)];

problem.grid.nTrajPts = 50;

end