function problem = quadrotor_1d_move_to(config)
addpath("..\models\")
%%% setup dynamics
model = quadrotor1d(config);
problem.model = model;
nState = model.dim.nState;
nControl = model.dim.nControl;

%%% setup cost
problem.cost.type = "LQR";
problem.cost.R = 1 * eye(nControl);
problem.cost.Q = 10 * eye(nState);

%%% setup constraints
% boundaries constraints
problem.constraints.boundarys.initState = [0;0];
problem.constraints.boundarys.finalState = [1;0];
problem.constraints.boundarys.finalTime = 1;

% bound constraints
problem.constraints.bounds.control.lower = -20*ones(nControl,1);
problem.constraints.bounds.control.upper = 20*ones(nControl,1);
problem.constraints.bounds.state.lower = -inf*ones(nState,1);
problem.constraints.bounds.state.upper = inf*ones(nState,1);
problem.grid.nTrajPts = 50;

end