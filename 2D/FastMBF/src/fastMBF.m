function Factors = fastMBF(N, fun, xx, xbox, kk, kbox, NG, tol)
    Nk = N;
    levels = floor(log2(Nk/NG));
    iter=1;
    Factors = cell(levels+1,2);
    kkidglobal = 1:size(kk,1);
    while(size(kk,1)>4*NG^2)
        ckbox = kbox/2;
        ck1s = ckbox(1,1);
        ck1e = ckbox(1,2);
        ck2s = ckbox(2,1);
        ck2e = ckbox(2,2);
        kkid = find( (kk(:,1)<ck1s | kk(:,1)>=ck1e) ...
            | (kk(:,2)<ck2s | kk(:,2)>=ck2e) );
        ckkid = find( kk(:,1)>=ck1s & kk(:,1)<ck1e ...
            & kk(:,2)>=ck2s & kk(:,2)<ck2e );
        Factors{iter,2} = kkidglobal(kkid);
        Factors{iter,1} = fastBF(fun, xx, xbox, kk(kkid,:), kbox, NG, tol);
        iter = iter+1;
        kbox = ckbox;
        kk = kk(ckkid,:);
        kkidglobal = kkidglobal(ckkid);
        Nk = Nk/2;
    end
    Factors{iter,2} = kkidglobal;
    Factors{iter,1} = fun(xx,kk);
end