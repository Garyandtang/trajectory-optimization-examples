function problem = cart_pole_swing_up(config)
addpath("..\models\")
%%% setup dynamics
model = cartPoleModel(config);
problem.model = model;
nState = model.dim.nState;
nControl = model.dim.nControl;

%%% setup cost
problem.cost.type = "LQR";
problem.cost.R = 1 * eye(nControl);
problem.cost.Q = 0 * eye(nState);

%%% setup constraints
% boundaries constraints
problem.constraints.boundarys.initState = zeros(nState,1);
problem.constraints.boundarys.finalState = [0;pi;0;0];
problem.constraints.boundarys.finalTime = 3;

% bound constraints
problem.constraints.bounds.control.lower = -inf*ones(nControl,1);
problem.constraints.bounds.control.upper = inf*ones(nControl,1);
problem.constraints.bounds.state.lower = -inf*ones(nState,1);
problem.constraints.bounds.state.upper = inf*ones(nState,1);
problem.grid.nTrajPts = 100;

end