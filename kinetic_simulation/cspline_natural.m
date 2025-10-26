function y2 = cspline_natural(y,dx)
%  给定数组y, 利用Compact_difference计算y的2阶导数 (其中边界处2阶导数设为0)
%  @param   n       数组长度
%  @param   y       数组
%  @param   dx      步长
%  @param   y2      返回2阶导数
% !>------------------------------------------------------------------------------------------------------
    n = length(y);
    y2 = zeros(n,1);
    y_prime = zeros(n,1);
    a = zeros(n,1);
    b = zeros(n,1);
    c = zeros(n,1);
    idx2m6 = 6.0/dx/dx;
    y_prime(2:n-1)  = (y(3:n)-2.0*y(2:n-1)+y(1:n-2))*idx2m6;
    y_prime(1) = 0.0;
    y_prime(n) = 0.0;
    a(1:n) = 1.0;
    c(1:n) = 1.0;
    b(1:n) = 4.0;
    y2 = tridiag(a,b,c,y_prime);
return;