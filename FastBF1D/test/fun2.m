function res = fun2(N,x,k)
    xk = x*k';
    tmp = (2*pi)* (xk);
    res = complex(cos(tmp),sin(tmp));
end