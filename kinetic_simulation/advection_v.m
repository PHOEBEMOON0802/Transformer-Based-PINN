function advection_v()
    global prm
    global v f ef g
    for i=2:prm.n+1
        for j=1:2*prm.m+1
            g(i,j) = splint_v(v(j)-prm.charge/prm.mass*ef(i)*prm.dt,i);
         end
    end
	for j = 1:prm.m+1
        g(prm.n+2,j)  = g(2,j);
		g(1,j) = g(prm.n+1,j);
    end
return;