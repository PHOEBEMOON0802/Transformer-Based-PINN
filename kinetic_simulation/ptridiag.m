function x=ptridiag(a, b, c, d)
    n = length(d);
    x = zeros(n,1);
    u = zeros(n,1);
    x2 = zeros(n,1);
    
    alpha = a(1);
    beta =  c(n);
    gamma = -b(1);
    a1 = a;
    b1 = b;
    c1 = c;
    b1(1) = b1(1) - gamma;
    b1(n) = b1(n) - alpha*beta/gamma;
    x = tridiag(a1,b1,c1,d);
    u = 0.d0 ; 
    u(1) = gamma ;
    u(n) = beta;
    b(1) = b(1) - gamma;
    b(n) = b(n) - alpha*beta/gamma;
    x2 = tridiag(a,b,c,u);
    factor = ( x(1)+ x(n)*alpha/gamma ) / (1.d0 + x2(1) + x2(n)*alpha/gamma);
    x = x - factor * x2;
return;