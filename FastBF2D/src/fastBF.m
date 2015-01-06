function Factor = fastBF(fun,xx,xbox,kk,kbox,NG,tol)

grid = fbf_grid(NG);

Nx = floor(sqrt(size(xx,1)));
Nk = floor(sqrt(size(kk,1)));

npx1 = 2^ceil(log2(sqrt(Nx)));
npx2 = 2^ceil(log2(sqrt(Nx)));
npk1 = 2^ceil(log2(sqrt(Nk)));
npk2 = 2^ceil(log2(sqrt(Nk)));
levels = floor(log2(Nx/npx1/NG));
LagrangeMatCell = cell(2,1);
LagrangeMatCell{1} = fbf_Lagrange(grid,grid/2);
LagrangeMatCell{2} = fbf_Lagrange(grid,grid/2+1/2);
TensorLagrangeMatCell = cell(2,2);
for a = 1:2
for b = 1:2
    TensorLagrangeMatCell{a,b} = kron(LagrangeMatCell{b},LagrangeMatCell{a});
end
end
Factor = struct('U',[],'GTol',[],'M',[],'HTol',[],'V',[]);

x1box = xbox(1,:);
x2box = xbox(2,:);
k1box = kbox(1,:);
k2box = kbox(2,:);

%---------------------------------------------------------------
%   Middle level construction

Mcell = cell(npx1,npx2,npk1,npk2);
LRU = cell(npx1,npx2,npk1,npk2);
LRV = cell(npx1,npx2,npk1,npk2);

for x1 = 1:npx1
for x2 = 1:npx2
    for k1 = 1:npk1
    for k2 = 1:npk2
        x1len = (x1box(2)-x1box(1))/npx1;
        x1s = x1box(1)+(x1-1)*x1len;
        x2len = (x2box(2)-x2box(1))/npx2;
        x2s = x2box(1)+(x2-1)*x2len;
        k1len = (k1box(2)-k1box(1))/npk1;
        k1s = k1box(1)+(k1-1)*k1len;
        k2len = (k2box(2)-k2box(1))/npk2;
        k2s = k2box(1)+(k2-1)*k2len;
        [x1grid,x2grid] = ndgrid(grid*x1len+x1s,grid*x2len+x2s);
        [k1grid,k2grid] = ndgrid(grid*k1len+k1s,grid*k2len+k2s);
        xgrid = [x1grid(:) x2grid(:)];
        kgrid = [k1grid(:) k2grid(:)];
        [Utmp,Stmp,Vtmp] = svdtrunc(fun(xgrid,kgrid),tol);
        Mcell{x1,x2,k1,k2} = size(Stmp,1);
        LRU{x1,x2,k1,k2} = Utmp*sqrt(Stmp);
        LRV{k1,k2,x1,x2} = sqrt(Stmp)*Vtmp';
    end
    end
end
end

offsetH = zeros(npx1,npx2,npk1,npk2);
sizeH = zeros(npx1,npx2,npk1,npk2);
for x1=1:npx1
for x2=1:npx2
    for k1=1:npk1
    for k2=1:npk2
        sizeH(x1,x2,k1,k2) = Mcell{x1,x2,k1,k2};
        if(k2>1)
            offsetH(x1,x2,k1,k2) = offsetH(x1,x2,k1,k2-1)+sizeH(x1,x2,k1,k2);
        elseif(k1>1)
            offsetH(x1,x2,k1,k2) = offsetH(x1,x2,k1-1,end)+sizeH(x1,x2,k1,k2);
        elseif(x2>1)
            offsetH(x1,x2,k1,k2) = offsetH(x1,x2-1,end,end)+sizeH(x1,x2,k1,k2);
        elseif(x1>1)
            offsetH(x1,x2,k1,k2) = offsetH(x1-1,end,end,end)+sizeH(x1,x2,k1,k2);
        end
    end
    end
end
end
midoffsetH = offsetH;
midsizeH = sizeH;

offsetW = zeros(npk1,npk2,npx1,npx2);
sizeW = zeros(npk1,npk2,npx1,npx2);
for k1=1:npk1
for k2=1:npk2
    for x1=1:npx1
    for x2=1:npx2
        sizeW(k1,k2,x1,x2) = Mcell{x1,x2,k1,k2};
        if(x2>1)
            offsetW(k1,k2,x1,x2) = offsetW(k1,k2,x1,x2-1)+sizeW(k1,k2,x1,x2);
        elseif(x1>1)
            offsetW(k1,k2,x1,x2) = offsetW(k1,k2,x1-1,end)+sizeW(k1,k2,x1,x2);
        elseif(k2>1)
            offsetW(k1,k2,x1,x2) = offsetW(k1,k2-1,end,end)+sizeW(k1,k2,x1,x2);
        elseif(k1>1)
            offsetW(k1,k2,x1,x2) = offsetW(k1-1,end,end,end)+sizeW(k1,k2,x1,x2);
        end
    end
    end
end
end
midoffsetW = offsetW;
midsizeW = sizeW;

totalel = 0;
for x1=1:npx1
for x2=1:npx2
    for k1=1:npk1
    for k2=1:npk2
        totalel = totalel + Mcell{x1,x2,k1,k2};
    end
    end
end
end

offset = 0;
XT = zeros(totalel,1);
YT = zeros(totalel,1);
ST = zeros(totalel,1);
for x1=1:npx1
for x2=1:npx2
    for k1=1:npk1
    for k2=1:npk2
        Mlen = Mcell{x1,x2,k1,k2};
        X = offsetH(x1,x2,k1,k2)+(1:sizeH(x1,x2,k1,k2));
        Y = offsetW(k1,k2,x1,x2)+(1:sizeW(k1,k2,x1,x2));
        idx = offset+(1:Mlen);
        XT(idx) = X(:);
        YT(idx) = Y(:);
        ST(idx) = ones(Mlen,1);
        if(~isempty(idx))
            offset = idx(end);
        end
    end
    end
end
end
MSpr = sparse(XT,YT,ST);

Factor.M = MSpr;
clear MSpr;

%---------------------------------------------------------------
%   Left factors construction
%           G
GTol = cell(levels,1);
npxx1 = npx1;
npxx2 = npx2;
npkk1 = npk1;
npkk2 = npk2;
offsetW = midoffsetH;
sizeW = midsizeH;

for ell = 1:levels
    
    npxx1 = npxx1*2;
    npxx2 = npxx2*2;
    npkk1 = npkk1/2;
    npkk2 = npkk2/2;
    Gcur = cell(npxx1,npxx2,npkk1,npkk2);
    LRUcur = cell(npxx1,npxx2,npkk1,npkk2);
    
    for x1 = 1:npxx1
    for x2 = 1:npxx2
        
        x1len = (x1box(2)-x1box(1))/npxx1;
        x1s = x1box(1)+(x1-1)*x1len;
        x1par = floor((x1-1)/2)+1;
        x1parlen = (x1box(2)-x1box(1))/npxx1*2;
        x1pars = x1box(1)+(x1par-1)*x1parlen;

        x2len = (x2box(2)-x2box(1))/npxx2;
        x2s = x2box(1)+(x2-1)*x2len;
        x2par = floor((x2-1)/2)+1;
        x2parlen = (x2box(2)-x2box(1))/npxx2*2;
        x2pars = x2box(1)+(x2par-1)*x2parlen;

        [x1grid,x2grid] = ndgrid(grid*x1len+x1s,grid*x2len+x2s);
        [x1pargrid,x2pargrid] = ndgrid(grid*x1parlen+x1pars,grid*x2parlen+x2pars);
        xgrid = [x1grid(:) x2grid(:)];
        xpargrid = [x1pargrid(:) x2pargrid(:)];
        
        LagrangeMat1 = fbf_Lagrange(unique(x1pargrid),unique(x1grid));
        LagrangeMat2 = fbf_Lagrange(unique(x2pargrid),unique(x2grid));
        TensorLagrangeMat = kron(LagrangeMat2,LagrangeMat1);
        
        for k1 = 1:npkk1
        for k2 = 1:npkk2
            
            GcurSub = cell(2,2);
            for a=1:2
            for b=1:2
                k1child = 2*k1-2+a;
                k1cen = k1box(1)+(k1child-1/2)*(k1box(2)-k1box(1))/npkk1/2;
                k2child = 2*k2-2+b;
                k2cen = k2box(1)+(k2child-1/2)*(k2box(2)-k2box(1))/npkk2/2;
                kcen = [k1cen k2cen];
                
                GcurSub{a,b} = diag(fun(xgrid,kcen))*(TensorLagrangeMat.'*(diag(1./fun(xpargrid,kcen))*LRU{x1par,x2par,k1child,k2child}));
            end
            end
            
            [Utmp,Stmp,Vtmp] = svdtrunc([GcurSub{1,1} GcurSub{1,2} GcurSub{2,1} GcurSub{2,2}],tol);
            Gcur{x1,x2,k1,k2} = Vtmp';
            LRUcur{x1,x2,k1,k2} = Utmp*Stmp;
        end
        end
    end
    end
    
    LRU = LRUcur;
    
    offsetH = zeros(npxx1,npxx2,npkk1,npkk2);
    sizeH = zeros(npxx1,npxx2,npkk1,npkk2);
    totalel = 0;
    for x1=1:npxx1
    for x2=1:npxx2
        for k1=1:npkk1
        for k2=1:npkk2
            totalel = totalel + numel(Gcur{x1,x2,k1,k2});
            sizeH(x1,x2,k1,k2) = size(Gcur{x1,x2,k1,k2},1);
            if(k2>1)
                offsetH(x1,x2,k1,k2) = offsetH(x1,x2,k1,k2-1)+sizeH(x1,x2,k1,k2);
            elseif(k1>1)
                offsetH(x1,x2,k1,k2) = offsetH(x1,x2,k1-1,end)+sizeH(x1,x2,k1,k2);
            elseif(x2>1)
                offsetH(x1,x2,k1,k2) = offsetH(x1,x2-1,end,end)+sizeH(x1,x2,k1,k2);
            elseif(x1>1)
                offsetH(x1,x2,k1,k2) = offsetH(x1-1,end,end,end)+sizeH(x1,x2,k1,k2);
            end
        end
        end
    end
    end
    
    offset = 0;
    XT = zeros(totalel,1);
    YT = zeros(totalel,1);
    ST = zeros(totalel,1);
    for x1=1:npxx1
    for x2=1:npxx2
        x1par = floor((x1-1)/2)+1;
        x2par = floor((x2-1)/2)+1;
        for k1=1:npkk1
        for k2=1:npkk2
            cursize = 0;
            for k1child = 2*k1-1:2*k1
            for k2child = 2*k2-1:2*k2
                GH = sizeH(x1,x2,k1,k2);
                GW = sizeW(x1par,x2par,k1child,k2child);
                [X,Y] = meshgrid(offsetH(x1,x2,k1,k2)+(1:GH),offsetW(x1par,x2par,k1child,k2child)+(1:GW));
                X = X';
                Y = Y';
                idx = offset+(1:GH*GW);
                XT(idx) = X(:);
                YT(idx) = Y(:);
                ST(idx) = Gcur{x1,x2,k1,k2}(:,cursize+(1:GW));
                cursize = cursize+GW;
                if(~isempty(idx))
                    offset = idx(end);
                end
            end
            end
        end
        end
    end
    end
    GTol{ell} = sparse(XT,YT,ST);
    offsetW = offsetH;
    sizeW = sizeH;
    
end

Factor.GTol = GTol;
clear GTol LRUcur Gcur;

%           U

Ucell = cell(npxx1,npxx2,npkk1,npkk2);

xxidx = fbf_subidx(xx,x1box,x2box,npxx1,npxx2);

for x1 = 1:npxx1
for x2 = 1:npxx2
    
    xxsub = xx(xxidx{x1,x2},:);
    xx1sub = xxsub(:,1);
    xx2sub = xxsub(:,2);
    
    x1len = (x1box(2)-x1box(1))/npxx1;
    x1s = x1box(1)+(x1-1)*x1len;
    x1grid = grid*x1len+x1s;
    LagrangeMat1 = fbf_Lagrange(x1grid,unique(xx1sub));
    
    x2len = (x2box(2)-x2box(1))/npxx2;
    x2s = x2box(1)+(x2-1)*x2len;
    x2grid = grid*x2len+x2s;
    LagrangeMat2 = fbf_Lagrange(x2grid,unique(xx2sub));
    
    [x1grid,x2grid] = ndgrid(x1grid,x2grid);
    xgrid = [x1grid(:) x2grid(:)];
    
    TensorLagrangeMat = kron(LagrangeMat2,LagrangeMat1);
    
    
    for k1 = 1:npkk1
    for k2 = 1:npkk2
        
        k1cen = k1box(1)+(k1-1/2)*(k1box(2)-k1box(1))/npkk1;
        k2cen = k2box(1)+(k2-1/2)*(k2box(2)-k2box(1))/npkk2;
        kcen = [k1cen k2cen];
        Ucell{x1,x2,k1,k2} = diag(fun(xxsub,kcen))*(TensorLagrangeMat.'*(diag(1./fun(xgrid,kcen))*LRU{x1,x2,k1,k2}));
    end
    end
end
end

totalel = 0;
for x1 = 1:npxx1
for x2 = 1:npxx2
    for k1 = 1:npkk1
    for k2 = 1:npkk2
        totalel = totalel + numel(Ucell{x1,x2,k1,k2});
    end
    end
end
end

offset = 0;
XT = zeros(totalel,1);
YT = zeros(totalel,1);
ST = zeros(totalel,1);
for x1=1:npxx1
for x2=1:npxx2
    for k1=1:npkk1
    for k2=1:npkk2
        UH = length(xxidx{x1,x2});
        UW = sizeW(x1,x2,k1,k2);
        [X,Y] = meshgrid(xxidx{x1,x2},offsetW(x1,x2,k1,k2)+(1:UW));
        X = X';
        Y = Y';
        idx = offset+(1:UH*UW);
        XT(idx) = X(:);
        YT(idx) = Y(:);
        ST(idx) = Ucell{x1,x2,k1,k2};
        if(~isempty(idx))
            offset = idx(end);
        end
    end
    end
end
end
USpr = sparse(XT,YT,ST);

Factor.U = USpr;
clear USpr LRU;
%---------------------------------------------------------------
%   Right factors construction
%            H
HTol = cell(levels,1);
npxx1 = npx1;
npxx2 = npx2;
npkk1 = npk1;
npkk2 = npk2;
offsetH = midoffsetW;
sizeH = midsizeW;

for ell = 1:levels
    
    npxx1 = npxx1/2;
    npxx2 = npxx2/2;
    npkk1 = npkk1*2;
    npkk2 = npkk2*2;
    Hcur = cell(npkk1,npkk2,npxx1,npxx2);
    LRVcur = cell(npkk1,npkk2,npxx1,npxx2);
    
    for k1 = 1:npkk1
    for k2 = 1:npkk2
        
        k1len = (k1box(2)-k1box(1))/npkk1;
        k1s = k1box(1)+(k1-1)*k1len;
        k1par = floor((k1-1)/2)+1;
        k1parlen = (k1box(2)-k1box(1))/npkk1*2;
        k1pars = k1box(1)+(k1par-1)*k1parlen;
        
        k2len = (k2box(2)-k2box(1))/npkk2;
        k2s = k2box(1)+(k2-1)*k2len;
        k2par = floor((k2-1)/2)+1;
        k2parlen = (k2box(2)-k2box(1))/npkk2*2;
        k2pars = k2box(1)+(k2par-1)*k2parlen;
        
        [k1grid,k2grid] = ndgrid(grid*k1len+k1s,grid*k2len+k2s);
        [k1pargrid,k2pargrid] = ndgrid(grid*k1parlen+k1pars,grid*k2parlen+k2pars);
        kgrid = [k1grid(:) k2grid(:)];
        kpargrid = [k1pargrid(:),k2pargrid(:)];
        
        LagrangeMat1 = fbf_Lagrange(unique(k1pargrid),unique(k1grid));
        LagrangeMat2 = fbf_Lagrange(unique(k2pargrid),unique(k2grid));
        TensorLagrangeMat = kron(LagrangeMat2,LagrangeMat1);
        
        for x1 = 1:npxx1
        for x2 = 1:npxx2
            
            HcurSub = cell(2,2);
            for a=1:2
            for b=1:2
                x1child = 2*x1-2+a;
                x1cen = x1box(1)+(x1child-1/2)*(x1box(2)-x1box(1))/npxx1/2;
                x2child = 2*x2-2+b;
                x2cen = x2box(1)+(x2child-1/2)*(x2box(2)-x2box(1))/npxx2/2;
                xcen = [x1cen x2cen];
                HcurSub{a,b} = LRV{k1par,k2par,x1child,x2child}*diag(1./fun(xcen,kpargrid))*TensorLagrangeMat*diag(fun(xcen,kgrid));
            end
            end
            
            [Utmp,Stmp,Vtmp] = svdtrunc([HcurSub{1,1};HcurSub{1,2};HcurSub{2,1};HcurSub{2,2}],tol);
            Hcur{k1,k2,x1,x2} = Utmp;
            LRVcur{k1,k2,x1,x2} = Stmp*Vtmp';
        end
        end
    end
    end
    
    LRV = LRVcur;
    
    offsetW = zeros(npkk1,npkk2,npxx1,npxx2);
    sizeW = zeros(npkk1,npkk2,npxx1,npxx2);
    totalel = 0;
    for k1=1:npkk1
    for k2=1:npkk2
        for x1=1:npxx1
        for x2=1:npxx2
            totalel = totalel + numel(Hcur{k1,k2,x1,x2});
            sizeW(k1,k2,x1,x2) = size(Hcur{k1,k2,x1,x2},2);
            if(x2>1)
                offsetW(k1,k2,x1,x2) = offsetW(k1,k2,x1,x2-1)+sizeW(k1,k2,x1,x2);
            elseif(x1>1)
                offsetW(k1,k2,x1,x2) = offsetW(k1,k2,x1-1,end)+sizeW(k1,k2,x1,x2);
            elseif(k2>1)
                offsetW(k1,k2,x1,x2) = offsetW(k1,k2-1,end,end)+sizeW(k1,k2,x1,x2);
            elseif(k1>1)
                offsetW(k1,k2,x1,x2) = offsetW(k1-1,end,end,end)+sizeW(k1,k2,x1,x2);
            end
        end
        end
    end
    end
    
    offset = 0;
    XT = zeros(totalel,1);
    YT = zeros(totalel,1);
    ST = zeros(totalel,1);
    for k1=1:npkk1
    for k2=1:npkk2
        k1par = floor((k1-1)/2)+1;
        k2par = floor((k2-1)/2)+1;
        for x1=1:npxx1
        for x2=1:npxx2
            cursize = 0;
            for x1child = 2*x1-1:2*x1
            for x2child = 2*x2-1:2*x2
                HH = sizeH(k1par,k2par,x1child,x2child);
                HW = sizeW(k1,k2,x1,x2);
                [X,Y] = meshgrid(offsetH(k1par,k2par,x1child,x2child)+(1:HH),offsetW(k1,k2,x1,x2)+(1:HW));
                X = X';
                Y = Y';
                idx = offset+(1:HH*HW);
                XT(idx) = X(:);
                YT(idx) = Y(:);
                ST(idx) = Hcur{k1,k2,x1,x2}(cursize+(1:HH),:);
                cursize = cursize+HH;
                if(~isempty(idx))
                    offset = idx(end);
                end
            end
            end
        end
        end
    end
    end
    HTol{ell} = sparse(XT,YT,ST);
    offsetH = offsetW;
    sizeH = sizeW;
    
end

Factor.HTol = HTol;
clear HTol LRVcur;

%           V

Vcell = cell(npkk1,npkk2,npxx1,npxx2);

kkidx = fbf_subidx(kk,k1box,k2box,npkk1,npkk2);

for k1 = 1:npkk1
for k2 = 1:npkk2
    
    kksub = kk(kkidx{k1,k2},:);
    kk1sub = kksub(:,1);
    kk2sub = kksub(:,2);
    
    k1len = (k1box(2)-k1box(1))/npkk1;
    k1s = k1box(1)+(k1-1)*k1len;
    k1grid = grid*k1len+k1s;
    LagrangeMat1 = fbf_Lagrange(k1grid,unique(kk1sub));
    
    k2len = (k2box(2)-k2box(1))/npkk2;
    k2s = k2box(1)+(k2-1)*k2len;
    k2grid = grid*k2len+k2s;
    LagrangeMat2 = fbf_Lagrange(k2grid,unique(kk2sub));
    
    [k1grid,k2grid] = ndgrid(k1grid,k2grid);
    kgrid = [k1grid(:) k2grid(:)];
    TensorLagrangeMat = kron(LagrangeMat2,LagrangeMat1);
    
    for x1 = 1:npxx1
    for x2 = 1:npxx2
        x1cen = x1box(1)+(x1-1/2)*(x1box(2)-x1box(1))/npxx1;
        x2cen = x2box(1)+(x2-1/2)*(x2box(2)-x2box(1))/npxx2;
        xcen = [x1cen x2cen];
        Vcell{k1,k2,x1,x2} = LRV{k1,k2,x1,x2}*diag(1./fun(xcen,kgrid))*TensorLagrangeMat*diag(fun(xcen,kksub));
    end
    end
end
end

totalel = 0;
for k1 = 1:npkk1
for k2 = 1:npkk2
    for x1 = 1:npxx1
    for x2 = 1:npxx2
        totalel = totalel + numel(Vcell{k1,k2,x1,x2});
    end
    end
end
end

offset = 0;
XT = zeros(totalel,1);
YT = zeros(totalel,1);
ST = zeros(totalel,1);
for k1=1:npkk1
for k2=1:npkk2
    for x1=1:npxx1
    for x2=1:npxx2
        VH = sizeH(k1,k2,x1,x2);
        VW = length(kkidx{k1,k2});
        [X,Y] = meshgrid(offsetH(k1,k2,x1,x2)+(1:VH),kkidx{k1,k2});
        X = X';
        Y = Y';
        idx = offset+(1:VH*VW);
        XT(idx) = X(:);
        YT(idx) = Y(:);
        ST(idx) = Vcell{k1,k2,x1,x2};
        if(~isempty(idx))
            offset = idx(end);
        end
    end
    end
end
end
VSpr = sparse(XT,YT,ST);

Factor.V = VSpr;
clear VSpr LRV;

end