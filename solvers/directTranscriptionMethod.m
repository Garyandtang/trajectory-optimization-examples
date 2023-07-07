function soln = directTranscriptionMethod(problem, config)
import casadi.*
addpath("../");

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
f1 = model.CTDynamic.evaluation.firstOrder; % this is 1st-order dynamics
f2 = model.CTDynamic.evaluation.secondOrder;
dt = finalTime / (nTrajPts - 1);
control.lower = problem.constraints.bounds.control.lower;
control.upper = problem.constraints.bounds.control.upper;
state.lower = problem.constraints.bounds.state.lower;
state.upper = problem.constraints.bounds.state.upper;

% object function transcription configuration
objConfig.dt = dt;
objConfig.Q = problem.cost.Q;
objConfig.R = problem.cost.R;
objConfig.method = config.method.objAppro;

% init solution
% todo: more reasonable init solution
w0 = zeros((nTrajPts)*(nControl+nState),1);


%%% Formulate the NLP
% define nlp container
nlpContainer.obj = 0;
nlpContainer.decVar.Z = [];
nlpContainer.decVar.U = [];
nlpContainer.bounds.lbz = [];
nlpContainer.bounds.ubz = [];
nlpContainer.bounds.lbu = [];
nlpContainer.bounds.ubu = [];
nlpContainer.constraints.g = [];
nlpContainer.constraints.lbg = [];
nlpContainer.constraints.ubg = [];


% decision variable on Z and U and var bound constraints
for i = 1 : nTrajPts
    zk = MX.sym(['X_' num2str(1)], nState);
    uk = MX.sym(['U_' num2str(1)], nControl);
    nlpContainer.decVar.Z = [nlpContainer.decVar.Z, zk];     
    nlpContainer.decVar.U = [nlpContainer.decVar.U, uk];    
    nlpContainer.bounds.lbu = [nlpContainer.bounds.lbu; control.lower];   % concatenated control lower bound
    nlpContainer.bounds.ubu = [nlpContainer.bounds.ubu; control.upper];   % concatenated control upper bound
    nlpContainer.bounds.lbz = [nlpContainer.bounds.lbz; state.lower];     % concatenated state lower bound 
    nlpContainer.bounds.ubz = [nlpContainer.bounds.ubz; state.upper];     % concatenated state upper bound
end

% defect constraint
nlpContainer = setDefectConstraints(nlpContainer, model, dt, config);

% objective value
nlpContainer = setObjFunction(nlpContainer, model, objConfig);

% other constraints
% start state constraint
nlpContainer.constraints.g = [nlpContainer.constraints.g; nlpContainer.decVar.Z(:, 1)];                     % add init state bounds to g
nlpContainer.constraints.lbg = [nlpContainer.constraints.lbg; initState];        % add lower bound for init state
nlpContainer.constraints.ubg = [nlpContainer.constraints.ubg; initState];        % add upper bound for init state
% end state constraint
nlpContainer.constraints.g = [nlpContainer.constraints.g; nlpContainer.decVar.Z(:, end)];                     % add init state bounds to g
nlpContainer.constraints.lbg = [nlpContainer.constraints.lbg; finalState];        % add lower bound for init state
nlpContainer.constraints.ubg = [nlpContainer.constraints.ubg; finalState];        % add upper bound for init state


% concatenate all decVar
J = nlpContainer.obj;
w = [reshape(nlpContainer.decVar.Z, numel(nlpContainer.decVar.Z),1); ...
     reshape(nlpContainer.decVar.U, numel(nlpContainer.decVar.U),1)];
lbw = [nlpContainer.bounds.lbz; nlpContainer.bounds.lbu];
ubw = [nlpContainer.bounds.ubz; nlpContainer.bounds.ubu];

g = nlpContainer.constraints.g;
lbg = nlpContainer.constraints.lbg;
ubg = nlpContainer.constraints.ubg;


tic 
% Create an NLP solver
prob = struct('f', J, 'x', w, 'g', g);
solver = nlpsol('solver', 'ipopt', prob);

% Solve the NLP
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);
soln.solverTime = toc;

% pack the result
dim.nState = [nState, nTrajPts];
dim.nControl = [nControl, nTrajPts];
[~, xSoln, uSoln] = unpackResult(w_opt,dim);
tSoln = linspace(0, finalTime, nTrajPts);
soln.tSoln = tSoln;
soln.qSoln = xSoln(1:nConfig,:);
soln.dqSoln = xSoln(nConfig+1:end,:);
soln.uSoln = uSoln;
% recover ddq result based on second-order system dyns
soln.ddqSoln = full(f2(soln.qSoln, soln.dqSoln, uSoln));

soln.interp.u = @(tt)( interp1(tSoln',uSoln',tt')' );
soln.interp.ddq = @(tt)( interp1(tSoln', soln.ddqSoln', tt')');
soln.interp.dq = @(tt)(bSpline2(tSoln, soln.dqSoln, soln.ddqSoln,tt));
soln.interp.q = @(tt)(bSpline3(tSoln, soln.qSoln, soln.dqSoln, soln.ddqSoln, tt));
soln.interp.x = @(t)([soln.interp.q(t);soln.interp.dq(t)]);

if (isfield(problem, "optimalSoln"))
   soln.interp.sysDymError = @(t)(soln.interp.q(t) - problem.optimalSoln.qSoln(t));
elseif (isfield(problem, "trueSoln"))
    trueSoln.interp.q = @(t)(bSpline3(problem.trueSoln.tSoln, problem.trueSoln.qSoln, problem.trueSoln.dqSoln, problem.trueSoln.ddqSoln,t));
%     trueSoln.interp.q = @(t)(bSpline2(problem.trueSoln.tSoln, problem.trueSoln.qSoln, problem.trueSoln.dqSoln, t));
%     trueSoln.interp.q = @(t)(interp1(problem.trueSoln.tSoln',problem.trueSoln.qSoln',t')');
    soln.interp.sysDymError = @(t)(soln.interp.q(t) - trueSoln.interp.q(t));
else 
    warning("trueSoln is not provided, use estimated system dynamics error")
    soln.interp.sysDymError = @(t)(full(f2(soln.interp.q(t),soln.interp.dq(t),soln.interp.u(t)))- ...
                               soln.interp.ddq(t));
end



% use romberg quadrature to estimate the absolute dynamic error
absSysDymError = @(t)(abs(soln.interp.sysDymError(t)));
nSegment = nTrajPts-1;
quadTol = 1e-12;   %Compute quadrature to this tolerance  
for i=1:nSegment
    soln.info.sysDymError(:,i) = rombergQuadrature(absSysDymError,tSoln([i,i+1]),quadTol);
end

end
%=============== helper function ======================
function nlpContainer = setDefectConstraints(nlpContainer, model, dt, config)
    nState = size(nlpContainer.decVar.Z, 1);
    nTrajs = size(nlpContainer.decVar.Z, 2);
    nConfig = nState / 2;
    f = model.CTDynamic.evaluation.firstOrder; % this is 1st-order dynamics
    f2 = model.CTDynamic.evaluation.secondOrder; % this is 2nd-order dynamics
    for i = 2 : nTrajs
        ukPrev = nlpContainer.decVar.U(:, i-1);
        zkPrev = nlpContainer.decVar.Z(:, i-1);
        zk = nlpContainer.decVar.Z(:, i);
     
    
        % calcuate defect constraints besed on trapezoidal rule
        switch config.method.dynamics

            case "first_order_euler"
                dzkPrev = f(zkPrev, ukPrev);
                g = zk - zkPrev - dt*dzkPrev;
            case "second_order_euler"
                % decode q and q_dot from state z
                qkPrev = zkPrev(1:nConfig, :);
                vkPrev = zkPrev(nConfig+1:end, :);
                qk = zk(1:nConfig, :);
                vk = zk(nConfig+1:end, :);
                % set defect constraints
                ddqkPrev = f2(qkPrev,vkPrev, ukPrev);
                g1 = vk - vkPrev - dt*ddqkPrev;
                g2 = qk - qkPrev - dt*vkPrev - 0.5*dt*dt*ddqkPrev;
                g = [g1; g2];
            case "first_order_rk4"
                k1 = dt*f(zkPrev, ukPrev);
                k2 = dt*f(zkPrev + 0.5*k1, ukPrev);
                k3 = dt*f(zkPrev + 0.5*k2, ukPrev);
                k4 = dt*f(zkPrev + k3, ukPrev);
                g = zk - zkPrev - (1/6)*(k1 + 2*k2 + 2*k3 + k4);
            case "second_order_rk4"
                % decode q and q_dot from state z
                qkPrev = zkPrev(1:nConfig, :);
                vkPrev = zkPrev(nConfig+1:end, :);
                qk = zk(1:nConfig, :);
                vk = zk(nConfig+1:end, :);
                % set defect constraints
                k1 = dt*f2(qkPrev, vkPrev, ukPrev);
                k2 = dt*f2(qkPrev+0.5*dt*vkPrev, vkPrev + 0.5*k1, ukPrev);
                k3 = dt*f2(qkPrev+0.5*dt*vkPrev, vkPrev + 0.5*k2, ukPrev);
                k4 = dt*f2(qkPrev+dt*vkPrev , vkPrev + k3, ukPrev);
                g1 = vk - vkPrev - (1/6)*(k1 + 2*k2 + 2*k3 + k4);
                g2 = qk - qkPrev - dt*vkPrev - (dt/30)*(6*k1+5*k2+3*k3+k4);
                g = [g1;g2];
            otherwise
                warning('Unexpected config.method.dynamics.')
        end
        nlpContainer.constraints.g = [nlpContainer.constraints.g; g];
        nlpContainer.constraints.lbg = [nlpContainer.constraints.lbg; zeros(nState, 1)];
        nlpContainer.constraints.ubg = [nlpContainer.constraints.ubg; zeros(nState,1)];
    end
end

function nlpContainer = setObjFunction(nlpContainer, model, objConfig)
    dt = objConfig.dt;
    R = objConfig.R;
    Q = objConfig.Q;
    nTrajPts = size(nlpContainer.decVar.Z, 2);

    for i = 2 : nTrajPts
        ukPrev = nlpContainer.decVar.U(:, i-1);
        uk = nlpContainer.decVar.U(:, i);
        xkPrev = nlpContainer.decVar.Z(:, i-1);
        xk = nlpContainer.decVar.Z(:,i);
        % calculate cost function by tapezoidal rule
        switch objConfig.method
            case "trapzoid_explict"
                nlpContainer.obj = nlpContainer.obj + dt*(ukPrev'*R*ukPrev + uk'*R*uk)/2;
                nlpContainer.obj = nlpContainer.obj + dt*(xkPrev'*Q*xkPrev + xk'*Q*xk)/2;
            case "trapzoid_implict"
                nlpContainer.obj = nlpContainer.obj + dt*(ukPrev'*R*ukPrev + uk'*R*ukPrev + uk'*R*uk)/3;
%                 nlpContainer.obj = nlpContainer.obj + dt*(xkPrev'*Q*xkPrev + xk'*Q*xkPrev + xk'*Q*xk)/3;
            otherwise
                error('Unexpected config.method.dynamics.')
        end
        
    end
end


