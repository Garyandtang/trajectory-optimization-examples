function config = setBlockMoveConfig()
%%% provide configuration for dynamic and problem
% Number of points for initialization:
config.grid.nTrajPts = 25;
config.finalTime = 1; % cannot change for analytical solution

% Bounds:
config.initState = zeros(2,1);
config.finalState = [1, 0]';

% Compute an initial guess at a trajectory:
config.guess = computeBMGuess(config);

% flag
config.flag.animationOn = 0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   help function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function guess = computeBMGuess(config)
% guess = computeGuess(config)
%
% This function computes an initial guess at the trajectory for block move
n = config.grid.nTrajPts;
tEnd = config.finalTime;

t = linspace(0, tEnd, n);

% assmue the block moves with constant acceleration.
acce = 2; % 2*d/t^2
guess.control = acce * ones(n,1);
q = acce * t;   % position
dq = 0.5* acce * (t.*t);    % velocity
guess.state = [q;dq];
guess.time = t;
end