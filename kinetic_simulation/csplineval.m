functionyi = csplineval(xi,x,S)
%CSPLINEVAL  Evaluate cubic spline polynomials
%   yi = CSPLINE(xi,S) returns the value, at the entries of xi, of the cubic 
%   piecewise polynomial f  with coefficient matrix S
%Input:
%       xi entries interested
%       x  breaks
%       S  coefficient matrix
%Output:
%       yiyi=f(xi) computed from cubic spline polynomials
l = length(x) - 1;
lx = length(xi);
xs = reshape(xi,1,lx);
if lx, [~,index] = histc(xs,[-inf,x(2:l),inf]);
else index = ones(1,lx);
end
% adjust for troubles, like evaluation sites that are NaN or +-inf
infxs = find(xs==inf);
if ~isempty(infxs), index(infxs) = N; end
nogoodxs = find(index==0);
if ~isempty(nogoodxs), xi(nogoodxs) = NaN; index(nogoodxs) = 1; end
% now go to local coordinates ...
xs = xs-x(index);
% ... and apply Heun algorithm for multiplication
yi = S(index,1);
for i=2:4
yi = xs(:).*yi + S(index,i);
end
end
