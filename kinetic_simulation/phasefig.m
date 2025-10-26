global v x f

    [V, X] = meshgrid(v, x);
    [C, h] = contourf(X, V, f, 50);
    h.LineStyle = 'none';