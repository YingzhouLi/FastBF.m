function relerr = fbf_check(N,fun,f,u,NC)
    u = reshape(u,[N N]);
    app = zeros(NC,1);
    ext = zeros(NC,1);
    k = -N/2:N/2-1;
    [k1,k2] = ndgrid(k);
    kk = [k1(:) k2(:)];
    pp = fbf_k2p(N,kk);
    for g=1:NC
        x1 = floor(rand(1)*N);
        x2 = floor(rand(1)*N);
        app(g) = u(x1+1,x2+1);
        ext(g) = sum(fun([x1 x2]/N,fbf_p2k(N,pp))*(f(:)));
    end
    err = app-ext;
    relerr = norm(err)/norm(ext);
end