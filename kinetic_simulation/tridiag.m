function x = tridiag(a,b,c,d)
    n = length(d);
    x = zeros(n,1);
    c(1) = c(1) / b(1);
    for i=2:n-1
       c(i) = c(i)/(b(i)-c(i-1)*a(i));
    end
    d(1) = d(1)/b(1);
    for i=2:n
       d(i) = (d(i)-d(i-1)*a(i))/(b(i)-c(i-1)*a(i));
    end
    x(n) = d(n);
    for i=n-1:-1:1
       x(i) = d(i) - c(i)*x(i+1);
    end
return;