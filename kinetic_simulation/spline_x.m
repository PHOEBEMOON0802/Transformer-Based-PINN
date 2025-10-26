function spline_x()
    global prm
    global f f2
    %!> 计算f(x,v)对x的2次导数
    %! the last column is set to be equal to the first
    %! to allow interpolation beyond xmin+Nx*dx (up to xmax!)
    f(prm.n+2,:) = f(2,:);
    f(1,:) = f(prm.n+1,:);
    for j = 1:2*prm.m+1
        f2(2:prm.n+1,j) = cspline_periodic(f(2:prm.n+1,j),prm.dx);
        f2(prm.n+2,j) = f2(2,j);
        f2(1,j) = f2(prm.n+1,j);
    end
return;