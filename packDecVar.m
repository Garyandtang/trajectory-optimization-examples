function [z,dim] = packDecVar(tSpan,x,u)
% [z,dim] = dimDecVar(t,x,u)
%
% This function collapses the time-span (t), state (x)
% and control (u) matricies into a single vector
%
% INPUTS:
%   tSpan = [1, 2] = [t0, tF] = time span
%   x = [ns, ms] = state matrix = [dimension in state space, grid point]
%   u = [nc, mc] = control matrix = [dimension in control space, grid point]
%
% OUTPUTS:
%   z = [2+ns*ms+nc*nc, 1] = vector of decision variables
%   dim = struct with the size of the state and control matricies
%       .nState = size(x);
%       .nControl = size(u);
%
% See Also: UNPACKDECVAR

dim.nState = size(x);
dim.nControl = size(u);

z = [tSpan(1);tSpan(2); reshape(x,numel(x),1); reshape(u,numel(u),1)];

end
