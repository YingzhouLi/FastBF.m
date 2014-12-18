function Factor = fastBF(fun,xx,xbox,kk,kbox,NG)

grid = fbf_grid(NG);

[Nx,~] = size(xx);
[Nk,~] = size(kk);

npx = 2^ceil(log2(sqrt(Nx))+0.5);
npk = 2^ceil(log2(sqrt(Nk))+0.5);
levels = floor(log2(Nx/npx/NG));
Factor = struct('U',[],'GTol',[],'M',[],'HTol',[],'V',[]);


%---------------------------------------------------------------
%   Middle level construction

Mcell = cell(npx,npk);

for x = 1:npx
    for k = 1:npk
        xlen = (xbox(2)-xbox(1))/npx;
        xs = xbox(1)+(x-1)*xlen;
        klen = (kbox(2)-kbox(1))/npk;
        ks = kbox(1)+(k-1)*klen;
        xgrid = grid*xlen+xs;
        kgrid = grid*klen+ks;
        Mcell{x,k} = fun(xgrid,kgrid);
    end
end

totalH = zeros(npx,1);
for x=1:npx
    for k=1:npk
        totalH(x) = totalH(x) + size(Mcell{x,k},1);
    end
end
currentH = zeros(npx,1);
for x=1:npx
    if(x>1)
        currentH(x) = currentH(x-1)+totalH(x-1);
    end
end

totalW = zeros(npk,1);
for k=1:npk
    for x=1:npx
        totalW(k) = totalW(k) + size(Mcell{x,k},2);
    end
end
currentW = zeros(npk,1);
for k=1:npk
    if(k>1)
        currentW(k) = currentW(k-1)+totalW(k-1);
    end
end

totalel = 0;
for x=1:npx
    for k=1:npk
        totalel = totalel + numel(Mcell{x,k});
    end
end

offset = 0;
XT = zeros(totalel,1);
YT = zeros(totalel,1);
ST = zeros(totalel,1);
for x=1:npx
    localH = currentH(x);
    for k=1:npk
        [MH,MW] = size(Mcell{x,k});
        [X,Y] = meshgrid(localH+(1:MH),currentW(k)+(1:MW));
        X = X';
        Y = Y';
        idx = offset+(1:MH*MW);
        XT(idx) = X(:);
        YT(idx) = Y(:);
        ST(idx) = Mcell{x,k};
        if(~isempty(idx))
            offset = idx(end);
        end
        localH = localH + MH;
        currentW(k) = currentW(k) + MW;
    end
end
MSpr = sparse(XT,YT,ST);

Factor.M = MSpr;
clear MSpr;

%---------------------------------------------------------------
%   Left factors construction
%           G

GTol = cell(levels,1);
npxx = npx;
npkk = npk;

for ell = 1:levels
    
    npxx = npxx*2;
    npkk = npkk/2;
    Gcur = cell(npxx,npkk);
    
    for x = 1:npxx
        
        xlen = (xbox(2)-xbox(1))/npxx;
        xs = xbox(1)+(x-1)*xlen;
        xpar = floor((x-1)/2)+1;
        xparlen = (xbox(2)-xbox(1))/npxx*2;
        xpars = xbox(1)+(xpar-1)*xparlen;
        xgrid = grid*xlen+xs;
        xpargrid = grid*xparlen+xpars;
        LagrangeMat = fbf_Lagrange(xpargrid,xgrid).';
        
        for k = 1:npkk
            kcen1 = kbox(1)+(2*k-3/2)*(kbox(2)-kbox(1))/npkk/2;
            GcurL = repmat(fun(xgrid,kcen1),[1 NG]).*LagrangeMat.*repmat((1./fun(xpargrid,kcen1)).',[NG 1]);
            
            kcen2 = kbox(1)+(2*k-1/2)*(kbox(2)-kbox(1))/npkk/2;
            GcurR = repmat(fun(xgrid,kcen2),[1 NG]).*LagrangeMat.*repmat((1./fun(xpargrid,kcen2)).',[NG 1]);
            
            Gcur{x,k} = [GcurL GcurR];
        end
    end
    
    totalel = 0;
    for x = 1:npxx
        for k = 1:npkk
            totalel = totalel + numel(Gcur{x,k});
        end
    end
    
    offset = 0;
    curH = 0;
    curW = 0;
    XT = zeros(totalel,1);
    YT = zeros(totalel,1);
    ST = zeros(totalel,1);
    for x=1:2:npxx
        localW = curW;
        for k=1:npkk
            [GH,GW] = size(Gcur{x,k});
            [X,Y] = meshgrid(curH+(1:GH),localW+(1:GW));
            X = X';
            Y = Y';
            idx = offset+(1:GH*GW);
            XT(idx) = X(:);
            YT(idx) = Y(:);
            ST(idx) = Gcur{x,k};
            if(~isempty(idx))
                offset = idx(end);
            end
            curH = curH + GH;
            localW = localW + GW;
        end
        localW = curW;
        for k=1:npkk
            [GH,GW] = size(Gcur{x+1,k});
            [X,Y] = meshgrid(curH+(1:GH),localW+(1:GW));
            X = X';
            Y = Y';
            idx = offset+(1:GH*GW);
            XT(idx) = X(:);
            YT(idx) = Y(:);
            ST(idx) = Gcur{x+1,k};
            if(~isempty(idx))
                offset = idx(end);
            end
            curH = curH + GH;
            localW = localW + GW;
        end
        curW = localW;
    end
    GTol{ell} = sparse(XT,YT,ST);
    
end

Factor.GTol = GTol;
clear GTol;

%           U

Ucell = cell(npxx,npkk);

for x = 1:npxx
    
    xidx = (1:Nx/npxx)+Nx/npxx*(x-1);
    xxsub = xx(xidx);
    xlen = (xbox(2)-xbox(1))/npxx;
    xs = xbox(1)+(x-1)*xlen;
    xgrid = grid*xlen+xs;
    LagrangeMat = fbf_Lagrange(xgrid,xxsub).';
    
    for k = 1:npkk
        kcen = kbox(1)+(k-1/2)*(kbox(2)-kbox(1))/npkk;
        Ucell{x,k} = repmat(fun(xxsub,kcen),[1 NG]).*LagrangeMat.*repmat((1./fun(xgrid,kcen)).',[size(xxsub,1) 1]);
    end
end

totalel = 0;
for x = 1:npxx
    for k = 1:npkk
        totalel = totalel + numel(Ucell{x,k});
    end
end

offset = 0;
curH = 0;
curW = 0;
XT = zeros(totalel,1);
YT = zeros(totalel,1);
ST = zeros(totalel,1);
for x=1:npxx
    for k=1:npkk
        [UH,UW] = size(Ucell{x,k});
        [X,Y] = meshgrid(curH+(1:UH),curW+(1:UW));
        X = X';
        Y = Y';
        idx = offset+(1:UH*UW);
        XT(idx) = X(:);
        YT(idx) = Y(:);
        ST(idx) = Ucell{x,k};
        if(~isempty(idx))
            offset = idx(end);
        end
        curW = curW + UW;
    end
    curH = curH + UH;
end
USpr = sparse(XT,YT,ST);

Factor.U = USpr;
clear USpr;
%---------------------------------------------------------------
%   Right factors construction
%            H

HTol = cell(levels,1);
npxx = npx;
npkk = npk;

for ell = 1:levels
    
    npxx = npxx/2;
    npkk = npkk*2;
    Hcur = cell(npkk,npxx);
    
    for k = 1:npkk
        
        klen = (kbox(2)-kbox(1))/npkk;
        ks = kbox(1)+(k-1)*klen;
        kpar = floor((k-1)/2)+1;
        kparlen = (kbox(2)-kbox(1))/npkk*2;
        kpars = kbox(1)+(kpar-1)*kparlen;
        kgrid = grid*klen+ks;
        kpargrid = grid*kparlen+kpars;
        LagrangeMat = fbf_Lagrange(kpargrid,kgrid);
        
        for x = 1:npxx
            xcen1 = xbox(1)+(2*x-3/2)*(xbox(2)-xbox(1))/npxx/2;
            HcurU = repmat((1./fun(xcen1,kpargrid)).',[1 NG]).*LagrangeMat.*repmat(fun(xcen1,kgrid),[NG 1]);
            
            xcen2 = xbox(1)+(2*x-1/2)*(xbox(2)-xbox(1))/npxx/2;
            HcurD = repmat((1./fun(xcen2,kpargrid)).',[1 NG]).*LagrangeMat.*repmat(fun(xcen2,kgrid),[NG 1]);
            
            Hcur{k,x} = [HcurU; HcurD];
        end
    end
    
    totalel = 0;
    for k = 1:npkk
        for x = 1:npxx
            totalel = totalel + numel(Hcur{k,x});
        end
    end
    
    offset = 0;
    curH = 0;
    curW = 0;
    XT = zeros(totalel,1);
    YT = zeros(totalel,1);
    ST = zeros(totalel,1);
    for k=1:2:npkk
        localH = curH;
        for x=1:npxx
            [HH,HW] = size(Hcur{k,x});
            [X,Y] = meshgrid(localH+(1:HH),curW+(1:HW));
            X = X';
            Y = Y';
            idx = offset+(1:HH*HW);
            XT(idx) = X(:);
            YT(idx) = Y(:);
            ST(idx) = Hcur{k,x};
            if(~isempty(idx))
                offset = idx(end);
            end
            localH = localH + HH;
            curW = curW + HW;
        end
        localH = curH;
        for x=1:npxx
            [HH,HW] = size(Hcur{k+1,x});
            [X,Y] = meshgrid(localH+(1:HH),curW+(1:HW));
            X = X';
            Y = Y';
            idx = offset+(1:HH*HW);
            XT(idx) = X(:);
            YT(idx) = Y(:);
            ST(idx) = Hcur{k+1,x};
            if(~isempty(idx))
                offset = idx(end);
            end
            localH = localH + HH;
            curW = curW + HW;
        end
        curH = localH;
    end
    HTol{ell} = sparse(XT,YT,ST);
    
end

Factor.HTol = HTol;
clear HTol;

%           V

Vcell = cell(npkk,npxx);

for k = 1:npkk
    
    kidx = (1:Nk/npkk)+Nk/npkk*(k-1);
    kksub = kk(kidx);
    klen = (kbox(2)-kbox(1))/npkk;
    ks = kbox(1)+(k-1)*klen;
    kgrid = grid*klen+ks;
    LagrangeMat = fbf_Lagrange(kgrid,kksub);
    
    for x = 1:npxx
        xcen = xbox(1)+(x-1/2)*(xbox(2)-xbox(1))/npxx;
        Vcell{k,x} = repmat((1./fun(xcen,kgrid)).',[1 size(kksub,1)]).*LagrangeMat.*repmat(fun(xcen,kksub),[NG 1]);
    end
end

totalel = 0;
for k = 1:npkk
    for x = 1:npxx
        totalel = totalel + numel(Vcell{k,x});
    end
end

offset = 0;
curH = 0;
curW = 0;
XT = zeros(totalel,1);
YT = zeros(totalel,1);
ST = zeros(totalel,1);
for k=1:npkk
    for x=1:npxx
        [VH,VW] = size(Vcell{k,x});
        [X,Y] = meshgrid(curH+(1:VH),curW+(1:VW));
        X = X';
        Y = Y';
        idx = offset+(1:VH*VW);
        XT(idx) = X(:);
        YT(idx) = Y(:);
        ST(idx) = Vcell{k,x};
        if(~isempty(idx))
            offset = idx(end);
        end
        curH = curH + VH;
    end
    curW = curW + VW;
end
VSpr = sparse(XT,YT,ST);

Factor.V = VSpr;
clear VSpr;

end