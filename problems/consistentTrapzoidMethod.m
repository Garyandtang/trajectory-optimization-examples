function soln = consistentTrapzoidMethod(problem)
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
f = model.CTDynamic.evaluation.secondOrder; % this is 2nd-order dynamics
dt = finalTime / (nTrajPts - 1);
control.lower = problem.constraints.bounds.control.lower;
control.upper = problem.constraints.bounds.control.upper;
state.lower = problem.constraints.bounds.state.lower;
state.upper = problem.constraints.bounds.state.upper;

w0 = ones((nTrajPts)*(nControl+nState),1);

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


% decision variable on Q V and U and var bound constraints
for i = 1 : nTrajPts
    qk = MX.sym(['Q_' num2str(i)], nConfig);
    vk = MX.sym(['V_' num2str(i)], nConfig);
    zk = [qk; vk];
    uk = MX.sym(['U_' num2str(i)], nControl);
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
    qkPrev = zkPrev(1:nConfig, :);
    vkPrev = zkPrev(nConfig+1:end, :);
    zk = nlpConstainer.decVar.Z(:, i);
    qk = zk(1:nConfig, :);
    vk = zk(nConfig+1:end, :);
    ddqkPrev = f(qkPrev,vkPrev, ukPrev);
    ddqk = f(qk, vk, uk);
    % calcuate defect constraints besed on trapezoidal rule
    g1 = vk - vkPrev -dt/2*(ddqk + ddqkPrev);
    g2 = qk - qkPrev -vkPrev*dt - dt^2*(2*ddqkPrev + ddqk)/6;
    nlpConstainer.constraints.g = [nlpConstainer.constraints.g; g1;g2];
    nlpConstainer.constraints.lbg = [nlpConstainer.constraints.lbg; zeros(nState, 1)];
    nlpConstainer.constraints.ubg = [nlpConstainer.constraints.ubg; zeros(nState,1)];
    % calculate cost function by tapezoidal rule
    nlpConstainer.obj = nlpConstainer.obj + dt*(ukPrev^2 + uk^2)/2;
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