function ret = splint_x(x_in,j)
    global prm
    global f f2
    
    xx = x_in;
    while(xx < 0 || xx >= prm.xl)
    	 if (xx < 0)
         	xx = xx + prm.xl;
         end
         if (xx>prm.xl)
        	xx = xx - prm.xl;
         end
    end
    ret = cspline2(f(2:prm.n+2,j), prm.dx, f2(2:prm.n+2,j), xx);
return;