function prm = initial
    global v x t xs
    global f f2 g
    global ef rho u p q t
    global fv 
    
    %rk0 = 0.01;
    rk0 = 0.5;
    alpha = 0.01;
    %alpha = 0.05;
    %alpha = 0.001;
    cte = 1.0/sqrt(2*pi);
    %vmax = 10;
    vmax = 3.0;
    xl = 2.0*pi/rk0;
    beta = 1;
    
    
    prm.mass = 1.0;
    prm.charge = -1.0;
    %prm.n = 64;
    %prm.m = 63;S
    prm.n = 200;
    prm.m = 64;
    prm.dt = 0.01;
    %prm.dt = 0.01;
    prm.nt = 2000;
    prm.xl = xl;
    prm.dv = vmax/prm.m;
    prm.dx = xl/prm.n;
    prm.vmax = vmax;
    v = (-prm.m:prm.m)*prm.dv;
    x = (-1:prm.n)*prm.dx;
    b = 2;
    a = sqrt(8);
    
    f = zeros(prm.n+2,2*prm.m+1);
    f2 = zeros(prm.n+2,2*prm.m+1);
    fv = zeros(2*prm.m+1,1);
    g = zeros(prm.n+2,2*prm.m+1);
    ef = zeros(prm.n+2,1);
    
    rho = zeros(prm.n+2,1);
    u = zeros(prm.n+2,1);
    p = zeros(prm.n+2,1);
    q = zeros(prm.n+2,1);
    t = zeros(prm.n+2,1); 
    xs = zeros(prm.n+2,1);
      
    for i = 2:prm.n+1
        pert = alpha*cos(rk0*x(i));
        for j = 1:2*prm.m+1
            f(i,j) = cte*exp(-v(j)*v(j)/2)*(1+pert);
        end
    end
    f(:,1) = 0;
    f(:,2*prm.m+1) = 0;
    f(prm.n+2,:) = f(2,:);
    f(1,:) = f(prm.n+1,:);
return