function Factor = fastBF(N,fun,xx,xbox,pp,pbox,NG,tol)

grid = fbf_grid(NG);

Nx = floor(sqrt(size(xx,1)));
Np = floor(sqrt(size(pp,1)));

npx1 = 2^ceil(log2(sqrt(Nx)));
npx2 = 2^ceil(log2(sqrt(Nx)));
npp1 = 2^ceil(log2(sqrt(Np)));
npp2 = 4*2^ceil(log2(sqrt(Np)));
levels = floor(log2(Nx/npx1/NG));
Factor = struct('U',[],'GTol',[],'M',[],'HTol',[],'V',[]);

x1box = xbox(1,:);
x2box = xbox(2,:);
p1box = pbox(1,:);
p2box = pbox(2,:);

%---------------------------------------------------------------
%   Middle level construction

Mcell = cell(npx1,npx2,npp1,npp2);
LRU = cell(npx1,npx2,npp1,npp2);
LRV = cell(npx1,npx2,npp1,npp2);

for x1 = 1:npx1
for x2 = 1:npx2
    for p1 = 1:npp1
    for p2 = 1:npp2
        x1len = (x1box(2)-x1box(1))/npx1;
        x1s = x1box(1)+(x1-1)*x1len;
        x2len = (x2box(2)-x2box(1))/npx2;
        x2s = x2box(1)+(x2-1)*x2len;
        p1len = (p1box(2)-p1box(1))/npp1;
        p1s = p1box(1)+(p1-1)*p1len;
        p2len = (p2box(2)-p2box(1))/npp2;
        p2s = p2box(1)+(p2-1)*p2len;
        [x1grid,x2grid] = ndgrid(grid*x1len+x1s,grid*x2len+x2s);
        [p1grid,p2grid] = ndgrid(grid*p1len+p1s,grid*p2len+p2s);
        xgrid = [x1grid(:) x2grid(:)];
        pgrid = [p1grid(:) p2grid(:)];
        [Utmp,Stmp,Vtmp] = svdtrunc(fun(xgrid,fbf_p2k(N,pgrid)),tol);
        Mcell{x1,x2,p1,p2} = size(Stmp,1);
        LRU{x1,x2,p1,p2} = Utmp*sqrt(Stmp);
        LRV{p1,p2,x1,x2} = sqrt(Stmp)*Vtmp';
    end
    end
end
end

offsetH = zeros(npx1,npx2,npp1,npp2);
sizeH = zeros(npx1,npx2,npp1,npp2);
for x1=1:npx1
for x2=1:npx2
    for p1=1:npp1
    for p2=1:npp2
        sizeH(x1,x2,p1,p2) = Mcell{x1,x2,p1,p2};
        if(p2>1)
            offsetH(x1,x2,p1,p2) = offsetH(x1,x2,p1,p2-1)+sizeH(x1,x2,p1,p2);
        elseif(p1>1)
            offsetH(x1,x2,p1,p2) = offsetH(x1,x2,p1-1,end)+sizeH(x1,x2,p1,p2);
        elseif(x2>1)
            offsetH(x1,x2,p1,p2) = offsetH(x1,x2-1,end,end)+sizeH(x1,x2,p1,p2);
        elseif(x1>1)
            offsetH(x1,x2,p1,p2) = offsetH(x1-1,end,end,end)+sizeH(x1,x2,p1,p2);
        end
    end
    end
end
end
midoffsetH = offsetH;
midsizeH = sizeH;

offsetW = zeros(npp1,npp2,npx1,npx2);
sizeW = zeros(npp1,npp2,npx1,npx2);
for p1=1:npp1
for p2=1:npp2
    for x1=1:npx1
    for x2=1:npx2
        sizeW(p1,p2,x1,x2) = Mcell{x1,x2,p1,p2};
        if(x2>1)
            offsetW(p1,p2,x1,x2) = offsetW(p1,p2,x1,x2-1)+sizeW(p1,p2,x1,x2);
        elseif(x1>1)
            offsetW(p1,p2,x1,x2) = offsetW(p1,p2,x1-1,end)+sizeW(p1,p2,x1,x2);
        elseif(p2>1)
            offsetW(p1,p2,x1,x2) = offsetW(p1,p2-1,end,end)+sizeW(p1,p2,x1,x2);
        elseif(p1>1)
            offsetW(p1,p2,x1,x2) = offsetW(p1-1,end,end,end)+sizeW(p1,p2,x1,x2);
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
    for p1=1:npp1
    for p2=1:npp2
        totalel = totalel + Mcell{x1,x2,p1,p2};
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
    for p1=1:npp1
    for p2=1:npp2
        Mlen = Mcell{x1,x2,p1,p2};
        X = offsetH(x1,x2,p1,p2)+(1:sizeH(x1,x2,p1,p2));
        Y = offsetW(p1,p2,x1,x2)+(1:sizeW(p1,p2,x1,x2));
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
nppp1 = npp1;
nppp2 = npp2;
offsetW = midoffsetH;
sizeW = midsizeH;

for ell = 1:levels
    
    npxx1 = npxx1*2;
    npxx2 = npxx2*2;
    nppp1 = nppp1/2;
    nppp2 = nppp2/2;
    Gcur = cell(npxx1,npxx2,nppp1,nppp2);
    LRUcur = cell(npxx1,npxx2,nppp1,nppp2);
    
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
        
        for p1 = 1:nppp1
        for p2 = 1:nppp2
            
            GcurSub = cell(2,2);
            for a=1:2
            for b=1:2
                p1child = 2*p1-2+a;
                p1cen = p1box(1)+(p1child-1/2)*(p1box(2)-p1box(1))/nppp1/2;
                p2child = 2*p2-2+b;
                p2cen = p2box(1)+(p2child-1/2)*(p2box(2)-p2box(1))/nppp2/2;
                pcen = [p1cen p2cen];
                
                GcurSub{a,b} = diag(fun(xgrid,fbf_p2k(N,pcen)))*(TensorLagrangeMat.'*(diag(1./fun(xpargrid,fbf_p2k(N,pcen)))*LRU{x1par,x2par,p1child,p2child}));
            end
            end
            
            [Utmp,Stmp,Vtmp] = svdtrunc([GcurSub{1,1} GcurSub{1,2} GcurSub{2,1} GcurSub{2,2}],tol);
            Gcur{x1,x2,p1,p2} = Vtmp';
            LRUcur{x1,x2,p1,p2} = Utmp*Stmp;
        end
        end
    end
    end
    
    LRU = LRUcur;
    
    offsetH = zeros(npxx1,npxx2,nppp1,nppp2);
    sizeH = zeros(npxx1,npxx2,nppp1,nppp2);
    totalel = 0;
    for x1=1:npxx1
    for x2=1:npxx2
        for p1=1:nppp1
        for p2=1:nppp2
            totalel = totalel + numel(Gcur{x1,x2,p1,p2});
            sizeH(x1,x2,p1,p2) = size(Gcur{x1,x2,p1,p2},1);
            if(p2>1)
                offsetH(x1,x2,p1,p2) = offsetH(x1,x2,p1,p2-1)+sizeH(x1,x2,p1,p2);
            elseif(p1>1)
                offsetH(x1,x2,p1,p2) = offsetH(x1,x2,p1-1,end)+sizeH(x1,x2,p1,p2);
            elseif(x2>1)
                offsetH(x1,x2,p1,p2) = offsetH(x1,x2-1,end,end)+sizeH(x1,x2,p1,p2);
            elseif(x1>1)
                offsetH(x1,x2,p1,p2) = offsetH(x1-1,end,end,end)+sizeH(x1,x2,p1,p2);
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
        for p1=1:nppp1
        for p2=1:nppp2
            cursize = 0;
            for p1child = 2*p1-1:2*p1
            for p2child = 2*p2-1:2*p2
                GH = sizeH(x1,x2,p1,p2);
                GW = sizeW(x1par,x2par,p1child,p2child);
                [X,Y] = meshgrid(offsetH(x1,x2,p1,p2)+(1:GH),offsetW(x1par,x2par,p1child,p2child)+(1:GW));
                X = X';
                Y = Y';
                idx = offset+(1:GH*GW);
                XT(idx) = X(:);
                YT(idx) = Y(:);
                ST(idx) = Gcur{x1,x2,p1,p2}(:,cursize+(1:GW));
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

Ucell = cell(npxx1,npxx2,nppp1,nppp2);

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
    
    
    for p1 = 1:nppp1
    for p2 = 1:nppp2
        
        p1cen = p1box(1)+(p1-1/2)*(p1box(2)-p1box(1))/nppp1;
        p2cen = p2box(1)+(p2-1/2)*(p2box(2)-p2box(1))/nppp2;
        pcen = [p1cen p2cen];
        Ucell{x1,x2,p1,p2} = diag(fun(xxsub,fbf_p2k(N,pcen)))*(TensorLagrangeMat.'*(diag(1./fun(xgrid,fbf_p2k(N,pcen)))*LRU{x1,x2,p1,p2}));
    end
    end
end
end

totalel = 0;
for x1 = 1:npxx1
for x2 = 1:npxx2
    for p1 = 1:nppp1
    for p2 = 1:nppp2
        totalel = totalel + numel(Ucell{x1,x2,p1,p2});
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
    for p1=1:nppp1
    for p2=1:nppp2
        UH = length(xxidx{x1,x2});
        UW = sizeW(x1,x2,p1,p2);
        [X,Y] = meshgrid(xxidx{x1,x2},offsetW(x1,x2,p1,p2)+(1:UW));
        X = X';
        Y = Y';
        idx = offset+(1:UH*UW);
        XT(idx) = X(:);
        YT(idx) = Y(:);
        ST(idx) = Ucell{x1,x2,p1,p2};
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
nppp1 = npp1;
nppp2 = npp2;
offsetH = midoffsetW;
sizeH = midsizeW;

for ell = 1:levels
    
    npxx1 = npxx1/2;
    npxx2 = npxx2/2;
    nppp1 = nppp1*2;
    nppp2 = nppp2*2;
    Hcur = cell(nppp1,nppp2,npxx1,npxx2);
    LRVcur = cell(nppp1,nppp2,npxx1,npxx2);
    
    for p1 = 1:nppp1
    for p2 = 1:nppp2
        
        p1len = (p1box(2)-p1box(1))/nppp1;
        p1s = p1box(1)+(p1-1)*p1len;
        p1par = floor((p1-1)/2)+1;
        p1parlen = (p1box(2)-p1box(1))/nppp1*2;
        p1pars = p1box(1)+(p1par-1)*p1parlen;
        
        p2len = (p2box(2)-p2box(1))/nppp2;
        p2s = p2box(1)+(p2-1)*p2len;
        p2par = floor((p2-1)/2)+1;
        p2parlen = (p2box(2)-p2box(1))/nppp2*2;
        p2pars = p2box(1)+(p2par-1)*p2parlen;
        
        [p1grid,p2grid] = ndgrid(grid*p1len+p1s,grid*p2len+p2s);
        [p1pargrid,p2pargrid] = ndgrid(grid*p1parlen+p1pars,grid*p2parlen+p2pars);
        pgrid = [p1grid(:) p2grid(:)];
        ppargrid = [p1pargrid(:),p2pargrid(:)];
        
        LagrangeMat1 = fbf_Lagrange(unique(p1pargrid),unique(p1grid));
        LagrangeMat2 = fbf_Lagrange(unique(p2pargrid),unique(p2grid));
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
                HcurSub{a,b} = LRV{p1par,p2par,x1child,x2child}*diag(1./fun(xcen,fbf_p2k(N,ppargrid)))*TensorLagrangeMat*diag(fun(xcen,fbf_p2k(N,pgrid)));
            end
            end
            
            [Utmp,Stmp,Vtmp] = svdtrunc([HcurSub{1,1};HcurSub{1,2};HcurSub{2,1};HcurSub{2,2}],tol);
            Hcur{p1,p2,x1,x2} = Utmp;
            LRVcur{p1,p2,x1,x2} = Stmp*Vtmp';
        end
        end
    end
    end
    
    LRV = LRVcur;
    
    offsetW = zeros(nppp1,nppp2,npxx1,npxx2);
    sizeW = zeros(nppp1,nppp2,npxx1,npxx2);
    totalel = 0;
    for p1=1:nppp1
    for p2=1:nppp2
        for x1=1:npxx1
        for x2=1:npxx2
            totalel = totalel + numel(Hcur{p1,p2,x1,x2});
            sizeW(p1,p2,x1,x2) = size(Hcur{p1,p2,x1,x2},2);
            if(x2>1)
                offsetW(p1,p2,x1,x2) = offsetW(p1,p2,x1,x2-1)+sizeW(p1,p2,x1,x2);
            elseif(x1>1)
                offsetW(p1,p2,x1,x2) = offsetW(p1,p2,x1-1,end)+sizeW(p1,p2,x1,x2);
            elseif(p2>1)
                offsetW(p1,p2,x1,x2) = offsetW(p1,p2-1,end,end)+sizeW(p1,p2,x1,x2);
            elseif(p1>1)
                offsetW(p1,p2,x1,x2) = offsetW(p1-1,end,end,end)+sizeW(p1,p2,x1,x2);
            end
        end
        end
    end
    end
    
    offset = 0;
    XT = zeros(totalel,1);
    YT = zeros(totalel,1);
    ST = zeros(totalel,1);
    for p1=1:nppp1
    for p2=1:nppp2
        p1par = floor((p1-1)/2)+1;
        p2par = floor((p2-1)/2)+1;
        for x1=1:npxx1
        for x2=1:npxx2
            cursize = 0;
            for x1child = 2*x1-1:2*x1
            for x2child = 2*x2-1:2*x2
                HH = sizeH(p1par,p2par,x1child,x2child);
                HW = sizeW(p1,p2,x1,x2);
                [X,Y] = meshgrid(offsetH(p1par,p2par,x1child,x2child)+(1:HH),offsetW(p1,p2,x1,x2)+(1:HW));
                X = X';
                Y = Y';
                idx = offset+(1:HH*HW);
                XT(idx) = X(:);
                YT(idx) = Y(:);
                ST(idx) = Hcur{p1,p2,x1,x2}(cursize+(1:HH),:);
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

Vcell = cell(nppp1,nppp2,npxx1,npxx2);

ppidx = fbf_subidx(pp,p1box,p2box,nppp1,nppp2);

for p1 = 1:nppp1
for p2 = 1:nppp2
    
    ppsub = pp(ppidx{p1,p2},:);
    pp1sub = ppsub(:,1);
    pp2sub = ppsub(:,2);
    
    p1len = (p1box(2)-p1box(1))/nppp1;
    p1s = p1box(1)+(p1-1)*p1len;
    p1grid = grid*p1len+p1s;
    LagrangeMat1 = fbf_Lagrange(p1grid,pp1sub);
    
    p2len = (p2box(2)-p2box(1))/nppp2;
    p2s = p2box(1)+(p2-1)*p2len;
    p2grid = grid*p2len+p2s;
    LagrangeMat2 = fbf_Lagrange(p2grid,pp2sub);
    
    [p1grid,p2grid] = ndgrid(p1grid,p2grid);
    pgrid = [p1grid(:) p2grid(:)];
    TensorLagrangeMat = kron(LagrangeMat2,ones(NG,1)).*kron(ones(NG,1),LagrangeMat1);
    
    for x1 = 1:npxx1
    for x2 = 1:npxx2
        x1cen = x1box(1)+(x1-1/2)*(x1box(2)-x1box(1))/npxx1;
        x2cen = x2box(1)+(x2-1/2)*(x2box(2)-x2box(1))/npxx2;
        xcen = [x1cen x2cen];
        Vcell{p1,p2,x1,x2} = LRV{p1,p2,x1,x2}*diag(1./fun(xcen,fbf_p2k(N,pgrid)))*TensorLagrangeMat*diag(fun(xcen,fbf_p2k(N,ppsub)));
    end
    end
end
end

totalel = 0;
for p1 = 1:nppp1
for p2 = 1:nppp2
    for x1 = 1:npxx1
    for x2 = 1:npxx2
        totalel = totalel + numel(Vcell{p1,p2,x1,x2});
    end
    end
end
end

offset = 0;
Hfix = 0;
XT = zeros(totalel,1);
YT = zeros(totalel,1);
ST = zeros(totalel,1);
for p1=1:nppp1
for p2=1:nppp2
    for x1=1:npxx1
    for x2=1:npxx2
        VH = sizeH(p1,p2,x1,x2);
        VW = length(ppidx{p1,p2});
        [X,Y] = meshgrid(offsetH(p1,p2,x1,x2)+(1:VH),ppidx{p1,p2});
        Hfix = max(Hfix, offsetH(p1,p2,x1,x2)+VH);
        X = X';
        Y = Y';
        idx = offset+(1:VH*VW);
        XT(idx) = X(:);
        YT(idx) = Y(:);
        ST(idx) = Vcell{p1,p2,x1,x2};
        if(~isempty(idx))
            offset = idx(end);
        end
    end
    end
end
end
VSpr = sparse(XT,YT,ST,Hfix,N^2);

Factor.V = VSpr;
clear VSpr LRV;

end