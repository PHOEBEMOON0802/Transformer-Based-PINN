function ret = splint_v(v_in,i)
    global prm
    global f f2
    if (v_in < -prm.vmax || v_in >= prm.vmax)
         ret = 0.d0;
         return;
    end
    ret= cspline2(f(i,:), prm.dv, f2(i,:), v_in+prm.vmax);
return;