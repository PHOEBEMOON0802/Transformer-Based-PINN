function efield3
    global prm 
    global f t time
    global ef rho v u p q rhoe
    d = zeros(prm.n,1);
    aa = zeros(prm.n,1);
    bb = zeros(prm.n,1);
    cc = zeros(prm.n,1);
    phi = zeros(prm.n,1);
    densm = 0.d0;
    densu = 0.d0;
    densp = 0.d0;
    densq = 0.d0;
    
    for i = 1:prm.n
        rho(i+1) = 0.0;
        u(i+1) = 0.0;
        p(i+1) = 0.0;
        q(i+1) = 0.0;
        

        for j = 1:2*prm.m+1
            rho(i+1) = rho(i+1)+f(i+1,j);
            u(i+1) = u(i+1) + (v(j)*f(i+1,j));
            
        end
        rho(i+1) = rho(i+1)*prm.dv;
        if rho(i+1) ~= 0
            u(i+1) = u(i+1) / rho(i+1)
        else
            u(i+1) = 0
        end
        u(i+1) = u(i+1)*prm.dv;
                
        for j = 1:2*prm.m+1
            p(i+1) = p(i+1) + (v(j)-u(i+1))^2*f(i+1,j);
            q(i+1) = q(i+1) + (v(j)-u(i+1))^3*f(i+1,j);
        end
        
        
        p(i+1) = p(i+1)*prm.dv;
        q(i+1) = q(i+1)*prm.dv;
        densm = densm + rho(i+1);
        %densu = densu + u(i+1);
        %densp = densp + p(i+1);
        %densq = densq + q(i+1);
    end 
    rhoe = rho;
    densm = densm/prm.n;
    for i = 1:prm.n
        rho(i+1) = prm.charge*rho(i+1)-prm.charge*densm;
        %u(i+1) = prm.charge*u(i+1)-prm.charge*densu;
        %p(i+1) = prm.charge*p(i+1)-prm.charge*densp;
        %q(i+1) = prm.charge*q(i+1)-prm.charge*densq;
    end
    dx2d12 = prm.dx^2/12.0;
    d(1)   = (rho(prm.n+1) + 10.0*rho(2) + rho(3))*dx2d12;
    d(prm.n) = (rho(prm.n) + 10.0*rho(prm.n+1)+rho(2))*dx2d12;
    for i = 1:prm.n
        d(i) = (rho(i) + 10.0*rho(i+1) + rho(i+2))*dx2d12; 
    end
    for i = 1:prm.n
        aa(i) = -1.0;
        bb(i) = 2.0;
        cc(i) = -1.0;
    end
    phi =  ptridiag(aa,bb,cc,d);
    tddx = 3.0/prm.dx;
    for i=2:prm.n-1
       d(i)=(phi(i+1)-phi(i-1))*tddx;
    end
    d(1)=(phi(2)-phi(prm.n))*tddx;
    d(prm.n)=(phi(1)-phi(prm.n-1))*tddx;
    for i = 1:prm.n
        aa(i) = -1.0;
        bb(i) = -4.0;
        cc(i) = -1.0;
    end
    U = ptridiag(aa,bb,cc,d);
    emm = 0.d0;
    for i = 1:prm.n
        emm = emm + U(i);
    end
    emm = emm/prm.n;
    for i = 1:prm.n
        ef(i+1) = U(i) - emm;
    end
    ef(prm.n+2) = ef(2);
    ef(1) = ef(prm.n+1);
return;
