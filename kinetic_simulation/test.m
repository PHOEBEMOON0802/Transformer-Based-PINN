    clear global
    warning off
    global prm f g t xs time x v
    

    global rho u p q ef
    
    prm = initial();
    time = 0.0;

    
    
    
    for it = 1:prm.nt
       
       time = time + prm.dt; 
       spline_x
       advection_x_semi
       f = g;
       efield3;
       
       spline_v
       advection_v
       f = g;
       
       spline_x
       advection_x_semi
       f = g;
       
       
    end

    figure(5);
    subplot(3,2,1);
    plot(xs(1:prm.n), rho(1:prm.n));title('n_x, t=12');
    subplot(3,2,2);
    plot(xs(1:prm.n), u(1:prm.n));title('u_x, t=12');
    subplot(3,2,3);
    plot(xs(1:prm.n), p(1:prm.n));title('p_x, t=12');
    subplot(3,2,4);
    plot(xs(1:prm.n), q(1:prm.n));title('q_x, t=12');
    subplot(3,2,5);
    plot(xs(1:prm.n), ef(1:prm.n));title('Ex_x, t=12');
