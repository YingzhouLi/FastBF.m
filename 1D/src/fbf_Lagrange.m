function res = fbf_Lagrange(gs,ts)
  NG = size(gs,1);
  NT = size(ts,1);

  tmp = zeros(NG,NT);
  for a=1:NG
      gud = [1:a-1 a+1:NG];
      for b=1:NT
          cur = (ts(b)-gs(gud))./(gs(a)-gs(gud));
          tmp(a,b) = prod(cur);
      end
  end
  res = tmp;
end
