function relerr = fbf_check(N,fun,f,u,NC)
  u = reshape(u,[N N]);
  app = zeros(NC,1);
  ext = zeros(NC,1);
  [k1,k2] = ndgrid((-N/2:N/2-1));
  ks = [k1(:) k2(:)];
  for g=1:NC
    x1 = floor(rand(1)*N);
    x2 = floor(rand(1)*N);
    app(g) = u(x1+1,x2+1);
    ext(g) = sum(fun([x1 x2]/N,ks)*(f(:)));
  end
  err = app-ext;
  relerr = norm(err)/norm(ext);
end