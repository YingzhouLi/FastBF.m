function res = funF(N,x,k)
    tmp = (2*pi)* (x*k');
    res = complex(cos(tmp),sin(tmp));
end