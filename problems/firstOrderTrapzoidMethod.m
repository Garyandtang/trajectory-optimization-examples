function soln = firstOrderTrapzoidMethod(problem)
import casadi.*
%%% first-order method
% get block move model
model = problem.model;
nState = model.dim.nState;
nControl = model.dim.nControl;
nConfig = model.dim.nConfig;

% get problem constraints
initState = problem.constraints.boundarys.initState;
finalState = problem.constraints.boundarys.finalState;
finalTime = problem.constraints.boundarys.finalTime;
nTrajPts = problem.grid.nTrajPts;
f = model.CTDynamic.evaluation.firstOrder; % this is 1st-order dynamics
dt = finalTime / (nTrajPts - 1);
control.lower = problem.constraints.bounds.control.lower;
control.upper = problem.constraints.bounds.control.upper;
state.lower = problem.constraints.bounds.state.lower;
state.upper = problem.constraints.bounds.state.upper;

w0 = zeros((nTrajPts)*(nControl+nState),1);

%%% Formulate the NLP
% define nlp container
nlpConstainer.obj = 0;
nlpConstainer.decVar.Z = [];
nlpConstainer.decVar.U = [];
nlpConstainer.bounds.lbz = [];
nlpConstainer.bounds.ubz = [];
nlpConstainer.bounds.lbu = [];
nlpConstainer.bounds.ubu = [];
nlpConstainer.constraints.g = [];
nlpConstainer.constraints.lbg = [];
nlpConstainer.constraints.ubg = [];


% decision variable on Z and U and var bound constraints
for i = 1 : nTrajPts
    zk = MX.sym(['X_' num2str(1)], nState);
    uk = MX.sym(['U_' num2str(1)], nControl);
    nlpConstainer.decVar.Z = [nlpConstainer.decVar.Z, zk];     
    nlpConstainer.decVar.U = [nlpConstainer.decVar.U, uk];    
    nlpConstainer.bounds.lbu = [nlpConstainer.bounds.lbu; control.lower];   % concatenated control lower bound
    nlpConstainer.bounds.ubu = [nlpConstainer.bounds.ubu; control.upper];   % concatenated control upper bound
    nlpConstainer.bounds.lbz = [nlpConstainer.bounds.lbz; state.lower];     % concatenated state lower bound 
    nlpConstainer.bounds.ubz = [nlpConstainer.bounds.ubz; state.upper];     % concatenated state upper bound
end

% defect constraints and objective value
for i = 2 : nTrajPts
    ukPrev = nlpConstainer.decVar.U(:, i-1);
    uk = nlpConstainer.decVar.U(:, i-1);
    zkPrev = nlpConstainer.decVar.Z(:, i-1);
    zk = nlpConstainer.decVar.Z(:, i);
    dzkPrev = f(zkPrev, ukPrev);
    dzk = f(zk, uk);
    % calcuate defect constraints besed on trapezoidal rule
    g1 = zk - zkPrev - dt/2*(dzk + zkPrev);
    nlpConstainer.constraints.g = [nlpConstainer.constraints.g; g1];
    nlpConstainer.constraints.lbg = [nlpConstainer.constraints.lbg; zeros(nState, 1)];
    nlpConstainer.constraints.ubg = [nlpConstainer.constraints.ubg; zeros(nState,1)];
    % calculate cost function by tapezoidal rule
    nlpConstainer.obj = nlpConstainer.obj + dt*(ukPrev^2 +ukPrev*uk+ uk^2)/3;
end

% other constraints
% start state constraint
nlpConstainer.constraints.g = [nlpConstainer.constraints.g; nlpConstainer.decVar.Z(:, 1)];                     % add init state bounds to g
nlpConstainer.constraints.lbg = [nlpConstainer.constraints.lbg; initState];        % add lower bound for init state
nlpConstainer.constraints.ubg = [nlpConstainer.constraints.ubg; initState];        % add upper bound for init state
% end state constraint
nlpConstainer.constraints.g = [nlpConstainer.constraints.g; nlpConstainer.decVar.Z(:, end)];                     % add init state bounds to g
nlpConstainer.constraints.lbg = [nlpConstainer.constraints.lbg; finalState];        % add lower bound for init state
nlpConstainer.constraints.ubg = [nlpConstainer.constraints.ubg; finalState];        % add upper bound for init state


% concatenate all decVar
J = nlpConstainer.obj;
w = [reshape(nlpConstainer.decVar.Z, numel(nlpConstainer.decVar.Z),1); ...
     reshape(nlpConstainer.decVar.U, numel(nlpConstainer.decVar.U),1)];
lbw = [nlpConstainer.bounds.lbz; nlpConstainer.bounds.lbu];
ubw = [nlpConstainer.bounds.ubz; nlpConstainer.bounds.ubu];

g = nlpConstainer.constraints.g;
lbg = nlpConstainer.constraints.lbg;
ubg = nlpConstainer.constraints.ubg;
% Create an NLP solver
prob = struct('f', J, 'x', w, 'g', g);
solver = nlpsol('solver', 'ipopt', prob);

% Solve the NLP
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);

soln = w_opt;
end