function soln = consistentTrapzoidMethod(problem, objApproximation)
import casadi.*
addpath("..")
addpath("../");
%%% first-order method
% get block move model
model = problem.model;
nState = model.dim.nState;
nControl = model.dim.nControl;
nConfig = model.dim.nConfig;
nTrajPts = problem.grid.nTrajPts;

% get problem constraints
initState = problem.constraints.boundarys.initState;
finalState = problem.constraints.boundarys.finalState;
finalTime = problem.constraints.boundarys.finalTime;

f = model.CTDynamic.evaluation.secondOrder; % this is 2nd-order dynamics
f1 = model.CTDynamic.evaluation.firstOrder;
dt = finalTime / (nTrajPts - 1);
control.lower = problem.constraints.bounds.control.lower;
control.upper = problem.constraints.bounds.control.upper;
state.lower = problem.constraints.bounds.state.lower;
state.upper = problem.constraints.bounds.state.upper;

w0 = ones((nTrajPts)*(nControl+nState),1);
tic
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
    uk = nlpConstainer.decVar.U(:, i);
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
    if (objApproximation == 1)
        nlpConstainer.obj = nlpConstainer.obj + dt*(ukPrev^2 + ukPrev*uk + uk^2)/3;
    else
        nlpConstainer.obj = nlpConstainer.obj + dt*(ukPrev^2 + uk^2)/2;
    end
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
toc
% pack the result
dim.nState = [nState, nTrajPts];
dim.nControl = [nControl, nTrajPts];
[~, xSoln, uSoln] = unpackResult(w_opt,dim);
tSoln = linspace(0, finalTime, nTrajPts);
soln.tSoln = tSoln;
soln.qSoln = xSoln(1:nConfig,:);
soln.dqSoln = xSoln(nConfig+1:end,:);
soln.uSoln = uSoln;

ddqSoln = full(f(soln.qSoln, soln.dqSoln, uSoln));
soln.interp.u = @(t)( interp1(tSoln',uSoln',t')' );
soln.interp.ddq = @(t)(interp1(tSoln', ddqSoln',t')');
soln.interp.dq = @(t)(bSpline2(tSoln, soln.dqSoln, ddqSoln,t));
soln.interp.q = @(t)(bSpline3(tSoln, soln.qSoln, soln.dqSoln, ddqSoln,t));
soln.interp.x = @(t)([soln.interp.q(t);soln.interp.dq(t)]);

soln.interp.sysDymError = @(t)(full(f1(soln.interp.x(t),soln.interp.u(t)))- ...
                               [soln.interp.dq(t);soln.interp.ddq(t)]);

% use romberg quadrature to estimate the absolute dynamic error
absSysDymError = @(t)(abs(soln.interp.sysDymError(t)));
nSegment = nTrajPts-1;
quadTol = 1e-12;   %Compute quadrature to this tolerance  
% soln.info.dynError = zeros(nConfig,nSegment);
for i=1:nSegment
    soln.info.sysDymError(:,i) = rombergQuadrature(absSysDymError,tSoln([i,i+1]),quadTol);
end
config.dyn.g = 9.81;
config.dyn.l = 1;
config.dyn.m1 = 1;
config.dyn.m2 = 0.1;
P.plotFunc = @(t,z)( drawCartPole(t,z,config.dyn) );
P.speed = 0.7;
P.figNum = 102;
t = linspace(soln.tSoln(1),soln.tSoln(end),250);
z = soln.interp.x(t);
if 0
   animate(t,z,P)
end

end