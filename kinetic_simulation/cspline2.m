function ret = cspline2(y,dx,y2,x)
%  给定数组y和它的2阶导数，通过插值返回x点处值.
%  @param   n       数组长度
%  @param   y       数组
%  @param   dx      步长
%  @param   y2      2阶导数
%  @param   x       所求数据点
%  @return  f(x)处的值
    i = floor(x/dx)+1;
    my_b = x/dx-(i-1);
    my_a = 1.d0-my_b;
    my_c = (my_a^3-my_a)*dx^2/6;
    my_d = (my_b^3-my_b)*dx^2/6;
    ret = y(i)*my_a + y(i+1)*my_b + y2(i)*my_c + y2(i+1)*my_d;
return;