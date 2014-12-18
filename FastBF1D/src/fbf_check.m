function relerr = fbf_check(N,fun,f,u,NC)
  app = zeros(NC,1);
  ext = zeros(NC,1);
  for g=1:NC
    x = floor(rand(1)*N+1);
    app(g) = u(x);
    ext(g) = sum(fun((x-1)/N,(-N/2:N/2-1)')*(f(:)));
  end
  err = app-ext;
  relerr = norm(err)/norm(ext);
end