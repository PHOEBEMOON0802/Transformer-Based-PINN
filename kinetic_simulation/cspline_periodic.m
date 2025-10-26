
function y2=cspline_periodic(y,dx)
%  给定数组y, 利用Compact_difference计算y的2阶导数 (周期边界条件)
%  @param   n       数组长度
%  @param   y       数组
% @param   dx      步长
%  @param   y2      返回2阶导数
% 
    n = length(y);
    y2 = zeros(n,1);
    y_prime = zeros(n,1);
    a = zeros(n,1);
    b = zeros(n,1);
    c = zeros(n,1);
    idx2m6 = 6/dx/dx;
    y_prime(2:n-1)  = (y(3:n)-2*y(2:n-1)+y(1:n-2))*idx2m6;
    y_prime(1) = (y(2)-2.0*y(1) + y(n))*idx2m6;
    y_prime(n) = (y(n-1)-2.0*y(n) + y(1))*idx2m6;
    a(1:n) = 1.0;
    c(1:n) = 1.0;
    b(1:n) = 4.0;
    y2=ptridiag(a,b,c,y_prime);
  return;