function [tSpan, x, u] = unpackResult(z, dim)
ns = dim.nState;  % ns = size(x)
nc = dim.nControl; % nc = size(u)

is = 1:(ns(1)*ns(2));
ic = 1:(nc(1)*nc(2));

x = reshape(z(is),ns(1),ns(2));
u = reshape(z(is(end)+ic),nc(1),nc(2));
tSpan = zeros(1,2);
end