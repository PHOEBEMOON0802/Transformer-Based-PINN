function advection_x_semi()
    global prm
    global x v f g
    
    
    for i=2:prm.n+1
        for j = 1:2*prm.m+1
           g(i,j) = splint_x(x(i)-prm.dt*v(j)*0.5,j);
        end
    end
	for j = 1:2*prm.m+1
        g(prm.n+2,j)  = g(2,j);
		g(1,j) = g(prm.n+1,j);
    end
return;