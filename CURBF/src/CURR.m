function [UR,sk] = CURR(fun,x,k,r)

rr = 5*r;
xlen = size(x,1);
klen = size(k,1);

sk = k(sort(randsample(klen,min(klen,r))),:);

px = x(randsample(xlen,min(xlen,rr)),:);

UR = fun(px,sk)\fun(px,k);

end