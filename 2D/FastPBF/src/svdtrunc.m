function [U,S,V] = svdtrunc(A,tol)
    [U,S,V] = svd(A,'econ');
    idx = find(diag(S)>tol*S(1,1));%*max(size(A)));
    U = U(:,idx);
    S = S(idx,idx);
    V = V(:,idx);
end