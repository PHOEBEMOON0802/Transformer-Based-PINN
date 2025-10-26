function spline_v()
    global prm
    global f f2
    for i=2:prm.n+2
        f2(i,:) = cspline_natural(f(i,:), prm.dv);
    end
return;