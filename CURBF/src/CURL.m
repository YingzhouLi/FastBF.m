function [CU,sx] = CURL(fun,x,k,r)

rr = 5*r;
xlen = size(x,1);
klen = size(k,1);

sx = x(sort(randsample(xlen,min(xlen,r))),:);

pk = k(randsample(klen,min(klen,rr)),:);

CU = fun(x,pk)/fun(sx,pk);

end