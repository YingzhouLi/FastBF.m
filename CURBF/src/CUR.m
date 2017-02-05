function [U,sx,sk] = CUR(fun,x,k,r)

rr = 5*r;
xlen = size(x,1);
klen = size(k,1);

sx = x(sort(randsample(xlen,min(xlen,r))),:);
sk = k(sort(randsample(klen,min(klen,r))),:);

px = [sx; x(randsample(xlen,min(xlen,rr)),:)];
pk = [sk; k(randsample(klen,min(klen,rr)),:)];

U = fun(px,sk)\fun(px,pk)/fun(sx,pk);

end