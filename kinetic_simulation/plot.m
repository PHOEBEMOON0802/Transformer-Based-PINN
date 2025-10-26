global v x f

    [V, X] = meshgrid(v, x);
    [C, h] = contourf(X, V, f, 100);
    h.LineStyle = 'none';