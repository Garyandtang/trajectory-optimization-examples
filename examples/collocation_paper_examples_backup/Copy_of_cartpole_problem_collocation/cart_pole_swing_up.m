function problem = cart_pole_swing_up(config)
addpath("..\models\")
%%% setup dynamics
model = cartPoleModel(config);
problem.model = model;
nState = model.dim.nState;
nControl = model.dim.nControl;

%%% setup cost
problem.cost.type = "LQR";
problem.cost.R = 5 * eye(nControl);
problem.cost.Q = 0 * eye(nState);

%%% setup constraints
% boundaries constraints
problem.constraints.boundarys.initState = zeros(nState,1);
problem.constraints.boundarys.finalState = [4;pi;0;0];
problem.constraints.boundarys.finalTime = 4;

% bound constraints
problem.constraints.bounds.control.lower = -20*ones(nControl,1);
problem.constraints.bounds.control.upper = 20*ones(nControl,1);
problem.constraints.bounds.state.lower = -inf*ones(nState,1);
problem.constraints.bounds.state.upper = inf*ones(nState,1);
problem.grid.nTrajPts = 40;

end