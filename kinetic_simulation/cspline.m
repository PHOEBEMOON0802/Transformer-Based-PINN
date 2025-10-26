function [yi, S] = cspline(x, y, xi, KC, dy0, dyN)
% CSPLINE finds the cubic splines for the input data points (x,y)
% Input: 
%       x = [x0 x1 ... xN], y = [y0 y1 ... yN], xi = interpolation points
%       KC = 1/2 for 1st/2nd derivatives on boundary specified
%       KC = 3 for 2nd derivative on boundary extrapolated
%       dy0 : initial derivative
%       dy1 : final derivative
% Output:
%       S(n,K); n = 1:N, k = 1,4 in descending order
ifnargin< 6, dyN = 0; end, if nargin< 5, dy0 = 0; end
ifnargin< 4, KC = 0; end
N = length(x) - 1;
A = zeros(N+1, N+1); b = zeros(N+1,1); S = zeros(N+1,4);
k = 1:N;
h(k) = x(k+1) - x(k);
dy(k) = (y(k+1)-y(k))./h(k);
% boundary condition
if KC <= 1 
    % 1st derivatives specified
A(1,1:2) = [2*h(1) h(1)]; b(1) = 3*(dy(1)-dy0);
A(N+1,N:N+1) = [h(N) 2*h(N)]; b(N+1) = 3*(dyN-dy(N));
elseif KC == 2
A(1,1) = 2; b(1)=dy0; A(N+1,N+1)=2; b(N+1)=dyN;
else
A(1,1:3) = [h(2) -h(1)-h(2) h(1)];
A(N+1,N-1:N+1) = [h(N) -h(N)-h(N-1) h(N-1)];
end
for m = 2:N
A(m,m-1:m+1) = [h(m-1) 2*(h(m-1)+h(m)) h(m)];
b(m) = 3 * (dy(m)-dy(m-1));
end
S(:,3) = A\b;
for m = 1:N
S(m,4) = (S(m+1,3)-S(m,3))/3.0/h(m);
   S(m,2) = dy(m)-h(m)/3.0*(S(m+1,3)+2*S(m,3));
S(m,1) = y(m);
end
S = S(1:N,4:-1:1);
yi = csplineval(xi,x,S);
end
