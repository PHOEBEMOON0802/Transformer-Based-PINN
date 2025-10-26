
x = load('x.txt');
t = load('t.txt');
n = load('n.txt');
u = load('u.txt');
p = load('p.txt');
q = load('q.txt');
E = load('E.txt');
f = load('f.txt');
data = horzcat(x, t, n, u, p, q, E);
   
save(['data_LD_6.mat'], 'data', '-v6');
%save(['f.mat'], 'f', '-v6');