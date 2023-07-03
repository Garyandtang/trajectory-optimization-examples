function problem = five_link_biped_walk(config)
addpath("..\models\")
%%% setup dynamics
model = fiveLinkBiped(config);
problem.model = model;
nState = model.dim.nState;
nControl = model.dim.nControl;

%%% setup cost
problem.cost.type = "LQR";
problem.cost.R = 1 * eye(nControl);
problem.cost.Q = 0 * eye(nState);

%%% setup constraints
% boundaries constraints
problem.constraints.boundarys.initState = [0.0848
    0.4941
   -0.1262
   -0.3469
   -0.3587
   -1.4975
   -0.6145
    0.0161
   -0.5520
   -2.0082
];
problem.constraints.boundarys.finalState = [-0.3587
   -0.3469
   -0.1262
    0.4941
    0.0848
   -2.1689
   -1.6314
    1.0293
   -1.4010
    3.0268];
problem.constraints.boundarys.finalTime = 0.8;

% bound constraints
problem.constraints.bounds.control.lower = -200*ones(nControl,1);
problem.constraints.bounds.control.upper = 200*ones(nControl,1);
problem.constraints.bounds.state.lower = -inf*ones(nState,1);
problem.constraints.bounds.state.upper = inf*ones(nState,1);

problem.grid.nTrajPts = 25;

end