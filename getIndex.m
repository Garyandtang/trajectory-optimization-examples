function [tIdx,xIdx,uIdx] = getIndex(dim)

% This function returns the indices to extract decision variables for use
% in linear constraints to be passed to fmincon.

nDecVar = 2 + dim.nState(1)*dim.nState(2) + ...
    dim.nControl(1)*dim.nControl(2);
[tIdx,xIdx,uIdx] = unPackDecVar(1:nDecVar,dim);

end