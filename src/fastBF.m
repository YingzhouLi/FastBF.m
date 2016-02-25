function Factor = fastBF(fun,xx,xbox,kk,kbox,NG,tol)

grid = Chey_grid(NG);

[Nx,Dim] = size(xx);
[Nk,~  ] = size(kk);
Nx = Nx^(1/Dim)*ones(1,Dim);
Nk = Nk^(1/Dim)*ones(1,Dim);

npx = 2.^ceil(log2(sqrt(Nx))+0.5);
npk = 2.^ceil(log2(sqrt(Nk))+0.5);
levels = max(0,min(floor(log2(Nx./npx./NG))));

xxboxidx = zeros(size(xx));
npxx = npx*2^levels;
for i = 1:Dim
    edges = linspace(xbox(1,i),xbox(2,i)+1/Nx(i),npxx(i)+1);
    [~,xxboxidx(:,i)] = histc(xx(:,i),edges);
end
[xxboxidx,xxidx] = sortrows(xxboxidx,Dim:-1:1);
[xxC,xxIA,~] = unique(xxboxidx,'rows','stable');
xxC = [xxC;zeros(1,size(xxC,2))];
xxIA = [xxIA;size(xx,1)+1];
for itx = 1:prod(npxx)
    x = idx2vec(npxx,itx);
    if any(x ~= xxC(itx,:))
        xxC = [xxC(1:itx-1,:);x;xxC(itx:end,:)];
        xxIA = [xxIA(1:itx-1);xxIA(itx);xxIA(itx:end)];
    end
end



kkboxidx = zeros(size(kk));
npkk = npk*2^levels;
for i = 1:Dim
    edges = linspace(kbox(1,i),kbox(2,i)+1/Nk(i),npkk(i)+1);
    [~,kkboxidx(:,i)] = histc(kk(:,i),edges);
end
[kkboxidx,kkidx] = sortrows(kkboxidx,Dim:-1:1);
[kkC,kkIA,~] = unique(kkboxidx,'rows','stable');
kkC = [kkC;zeros(1,size(kkC,2))];
kkIA = [kkIA;size(kk,1)+1];
for itk = 1:prod(npkk)
    k = idx2vec(npkk,itk);
    if any(k ~= kkC(itk,:))
        kkC = [kkC(1:itk-1,:);k;kkC(itk:end,:)];
        kkIA = [kkIA(1:itk-1);kkIA(itk);kkIA(itk:end)];
    end
end

%===============================================================
% The butterfly factorization is of the form
%  K = U Gl ... G2 G1 M H1 H2 ... Hl V
% where U is stored in 'U', G1 G2 ... Gl are stored in 'GTol',
% M is stored in 'M', H1 H2 ... Hl are stored in 'HTol',
% and V is stored in 'V'.
%===============================================================
Factor = struct('U',[],'GTol',[],'M',[],'HTol',[],'V',[]);

%---------------------------------------------------------------
%   M construction

Mcell = cell(prod(npx),prod(npk));

for itx = 1:prod(npx)
    x = idx2vec(npx,itx);
    for itk = 1:prod(npk)
        k = idx2vec(npk,itk);
        xgrid = fbf_grid(x,npx,grid,xbox);
        kgrid = fbf_grid(k,npk,grid,kbox);
        Mcell{itx,itk} = fun(xgrid,kgrid);
    end
end

%---------------------------------------------------------------
%   G construction

GTolcell = cell(levels,1);
for ell = 1:levels
    npxx = npx*2^ell;
    npkk_child = npk/2^(ell-1);
    GTolcell{ell} = cell(prod(npxx),prod(npkk_child));
    for itx = 1:prod(npxx)
        x = idx2vec(npxx,itx);
        [xgrid,~] = fbf_grid(x,npxx,grid,xbox);
        [xpargrid,xparLgrid] = fbf_grid(floor((x-1)/2)+1,npxx/2,grid,xbox);
        LagrangeMat = fbf_Lagrange(xparLgrid,xgrid).';
        for itk_child = 1:prod(npkk_child)
            k_child = idx2vec(npkk_child,itk_child);
            kcen_child = kbox(1,:)+(k_child-1/2).*(kbox(2,:)-kbox(1,:))./npkk_child;
            GTolcell{ell}{itx,itk_child} = ...
                diag(fun(xgrid,kcen_child))*(LagrangeMat*...
                (diag(1./fun(xpargrid,kcen_child))));
        end
    end
end

%---------------------------------------------------------------
%   U construction

npxx = npx*2^levels;
npkk = npk/2^levels;
Ucell = cell(prod(npxx),prod(npkk));

for itx = 1:prod(npxx)
    x = idx2vec(npxx,itx);
    [xgrid,xLgrid] = fbf_grid(x,npxx,grid,xbox);
    xxsub = xx(xxidx(xxIA(itx):xxIA(itx+1)-1),:);
    LagrangeMat = fbf_Lagrange(xLgrid,xxsub).';
    for itk = 1:prod(npkk)
        k = idx2vec(npkk,itk);
        kcen = kbox(1,:)+(k-1/2).*(kbox(2,:)-kbox(1,:))./npkk;
        Ucell{itx,itk} = ...
            diag(fun(xxsub,kcen))*(LagrangeMat*(diag(1./fun(xgrid,kcen))));
    end
end

%---------------------------------------------------------------
%   H construction

HTolcell = cell(levels,1);
for ell = 1:levels
    npkk = npk*2^ell;
    npxx_child = npx/2^(ell-1);
    HTolcell{ell} = cell(prod(npkk),prod(npxx_child));
    for itk = 1:prod(npkk)
        k = idx2vec(npkk,itk);
        kgrid = fbf_grid(k,npkk,grid,kbox);
        [kpargrid,kparLgrid] = fbf_grid(floor((k-1)/2)+1,npkk/2,grid,kbox);
        LagrangeMat = fbf_Lagrange(kparLgrid,kgrid);
        for itx_child = 1:prod(npxx_child)
            x_child = idx2vec(npxx_child,itx_child);
            xcen_child = xbox(1,:)+(x_child-1/2).*(xbox(2,:)-xbox(1,:))./npxx_child;
            HTolcell{ell}{itk,itx_child} = ...
                diag(1./fun(xcen_child,kpargrid))*LagrangeMat*...
                diag(fun(xcen_child,kgrid));
        end
    end
end

%---------------------------------------------------------------
%   V construction

npkk = npk*2^levels;
npxx = npx/2^levels;
Vcell = cell(prod(npkk),prod(npxx));

for itk = 1:prod(npkk)
    k = idx2vec(npkk,itk);
    [kgrid,kLgrid] = fbf_grid(k,npkk,grid,kbox);
    kksub = kk(kkidx(kkIA(itk):kkIA(itk+1)-1),:);
    LagrangeMat = fbf_Lagrange(kLgrid,kksub);

    for itx = 1:prod(npxx)
        x = idx2vec(npxx,itx);
        xcen = xbox(1,:)+(x-1/2).*(xbox(2,:)-xbox(1,:))./npxx;
        Vcell{itk,itx} = ...
            diag(1./fun(xcen,kgrid))*LagrangeMat*diag(fun(xcen,kksub));
    end
end

%==============================================================
% Fast butterfly factorization compression
% The compression is split into two parts: outward compression
% and inward compression

[Mcell,GTolcell,Ucell,HTolcell,Vcell] = fastBF_outComp(levels,npx,npk,tol,...
                                            Mcell,GTolcell,Ucell,HTolcell,Vcell);
%[Mcell,GTolcell,Ucell,HTolcell,Vcell] = fastBF_inComp(levels,npx,npk,tol,...
%                                            Mcell,GTolcell,Ucell,HTolcell,Vcell);

%==============================================================
% Sparse matrix assembling
%  K = U Gl ... G2 G1 M H1 H2 ... Hl V
% where U, Gi, M, Hi, and V are sparse matrices

%--------------------------------------------------------------
%   M assembling

totalel = 0;
for itx = 1:prod(npx)
    for itk = 1:prod(npk)
        totalel = totalel + numel(Mcell{itx,itk});
    end
end

totalH = 0;
offsetH = zeros(prod(npx),prod(npk));
for itx = 1:prod(npx)
    for itk = 1:prod(npk)
        offsetH(itx,itk) = totalH;
        totalH = totalH + size(Mcell{itx,itk},1);
    end
end

totalW = 0;
offsetW = zeros(prod(npk),prod(npx));
for itk = 1:prod(npk)
    for itx = 1:prod(npx)
        offsetW(itk,itx) = totalW;
        totalW = totalW + size(Mcell{itx,itk},2);
    end
end

XT = zeros(totalel,1);
YT = zeros(totalel,1);
ST = zeros(totalel,1);
offset = 0;
for itx = 1:prod(npx)
    for itk = 1:prod(npk)
        Mat = Mcell{itx,itk};
        [MatH,MatW] = size(Mat);
        [X,Y] = meshgrid(offsetH(itx,itk)+(1:MatH),offsetW(itk,itx)+(1:MatW));
        X = X';
        Y = Y';
        idx = offset+(1:MatH*MatW);
        XT(idx) = X(:);
        YT(idx) = Y(:);
        ST(idx) = Mat(:);
        if(~isempty(idx))
            offset = idx(end);
        end
    end
end
Factor.M = sparse(XT,YT,ST,totalH,totalW);


%---------------------------------------------------------------
%   G assembling

Factor.GTol = cell(levels,1);
for ell = 1:levels
    npxx = npx*2^ell;
    npkk = npk/2^ell;
    npxx_par = npx*2^(ell-1);
    npkk_child = npk/2^(ell-1);
    totalel = 0;
    for itx = 1:prod(npxx)
        for itk_child = 1:prod(npkk_child)
            totalel = totalel + numel(GTolcell{ell}{itx,itk_child});
        end
    end

    totalH = 0;
    offsetH = zeros(prod(npxx),prod(npkk_child));
    for itx = 1:prod(npxx)
        for itk = 1:prod(npkk)
            for it_child = 1:2^Dim
                itk_child = vec2idx(npkk_child,(idx2vec(npkk,itk)-1)*2+idx2vec(2*ones(1,Dim),it_child));
                offsetH(itx,itk_child) = totalH;
            end
            totalH = totalH + size(GTolcell{ell}{itx,itk_child},1);
        end
    end

    totalW = 0;
    offsetW = zeros(prod(npxx),prod(npkk_child));
    for itx_par = 1:prod(npxx_par)
        for itk_child = 1:prod(npkk_child)
            for it_child = 1:2^Dim
                itx = vec2idx(npxx,(idx2vec(npxx_par,itx_par)-1)*2+idx2vec(2*ones(1,Dim),it_child));
                offsetW(itx,itk_child) = totalW;
            end
            totalW = totalW + size(GTolcell{ell}{itx,itk_child},2);
        end
    end

    XT = zeros(totalel,1);
    YT = zeros(totalel,1);
    ST = zeros(totalel,1);
    offset = 0;
    for itx = 1:prod(npxx)
        for itk_child = 1:prod(npkk_child)
            Mat = GTolcell{ell}{itx,itk_child};
            [MatH,MatW] = size(Mat);
            [X,Y] = meshgrid(offsetH(itx,itk_child)+(1:MatH),offsetW(itx,itk_child)+(1:MatW));
            X = X';
            Y = Y';
            idx = offset+(1:MatH*MatW);
            XT(idx) = X(:);
            YT(idx) = Y(:);
            ST(idx) = Mat(:);
            if(~isempty(idx))
                offset = idx(end);
            end
        end
    end
    Factor.GTol{ell} = sparse(XT,YT,ST,totalH,totalW);
end

%---------------------------------------------------------------
%   U assembling

npxx = npx*2^levels;
npkk = npk/2^levels;
totalel = 0;
for itx = 1:prod(npxx)
    for itk = 1:prod(npkk)
        totalel = totalel + numel(Ucell{itx,itk});
    end
end

totalH = size(xx,1);

totalW = 0;
offsetW = zeros(prod(npxx),prod(npkk));
for itx = 1:prod(npxx)
    for itk = 1:prod(npkk)
        offsetW(itx,itk) = totalW;
        totalW = totalW + size(Ucell{itx,itk},2);
    end
end

XT = zeros(totalel,1);
YT = zeros(totalel,1);
ST = zeros(totalel,1);
offset = 0;
for itx = 1:prod(npxx)
    for itk = 1:prod(npkk)
        Mat = Ucell{itx,itk};
        [MatH,MatW] = size(Mat);
        [X,Y] = meshgrid(xxidx(xxIA(itx):xxIA(itx+1)-1),offsetW(itx,itk)+(1:MatW));
        X = X';
        Y = Y';
        idx = offset+(1:MatH*MatW);
        XT(idx) = X(:);
        YT(idx) = Y(:);
        ST(idx) = Mat(:);
        if(~isempty(idx))
            offset = idx(end);
        end
    end
end
Factor.U = sparse(XT,YT,ST,totalH,totalW);

%---------------------------------------------------------------
%   H assembling

Factor.HTol = cell(levels,1);
for ell = 1:levels
    npkk = npk*2^ell;
    npxx = npx/2^ell;
    npkk_par = npk*2^(ell-1);
    npxx_child = npx/2^(ell-1);
    totalel = 0;
    for itk = 1:prod(npkk)
        for itx_child = 1:prod(npxx_child)
            totalel = totalel + numel(HTolcell{ell}{itk,itx_child});
        end
    end

    totalW = 0;
    offsetW = zeros(prod(npkk),prod(npxx_child));
    for itk = 1:prod(npkk)
        for itx = 1:prod(npxx)
            for it_child = 1:2^Dim
                itx_child = vec2idx(npxx_child,(idx2vec(npxx,itx)-1)*2+idx2vec(2*ones(1,Dim),it_child));
                offsetW(itk,itx_child) = totalW;
            end
            totalW = totalW + size(HTolcell{ell}{itk,itx_child},2);
        end
    end

    totalH = 0;
    offsetH = zeros(prod(npkk),prod(npxx_child));
    for itk_par = 1:prod(npkk_par)
        for itx_child = 1:prod(npxx_child)
            for it_child = 1:2^Dim
                itk = vec2idx(npkk,(idx2vec(npkk_par,itk_par)-1)*2+idx2vec(2*ones(1,Dim),it_child));
                offsetH(itk,itx_child) = totalH;
            end
            totalH = totalH + size(HTolcell{ell}{itk,itx_child},1);
        end
    end

    XT = zeros(totalel,1);
    YT = zeros(totalel,1);
    ST = zeros(totalel,1);
    offset = 0;
    for itk = 1:prod(npkk)
        for itx_child = 1:prod(npxx_child)
            Mat = HTolcell{ell}{itk,itx_child};
            [MatH,MatW] = size(Mat);
            [X,Y] = meshgrid(offsetH(itk,itx_child)+(1:MatH),offsetW(itk,itx_child)+(1:MatW));
            X = X';
            Y = Y';
            idx = offset+(1:MatH*MatW);
            XT(idx) = X(:);
            YT(idx) = Y(:);
            ST(idx) = Mat(:);
            if(~isempty(idx))
                offset = idx(end);
            end
        end
    end
    Factor.HTol{ell} = sparse(XT,YT,ST,totalH,totalW);
end

%---------------------------------------------------------------
%   V assembling

npkk = npk*2^levels;
npxx = npx/2^levels;
totalel = 0;
for itk = 1:prod(npkk)
    for itx = 1:prod(npxx)
        totalel = totalel + numel(Vcell{itk,itx});
    end
end

totalH = 0;
offsetH = zeros(prod(npkk),prod(npxx));
for itk = 1:prod(npkk)
    for itx = 1:prod(npxx)
        offsetH(itk,itx) = totalH;
        totalH = totalH + size(Vcell{itk,itx},1);
    end
end

totalW = size(kk,1);

XT = zeros(totalel,1);
YT = zeros(totalel,1);
ST = zeros(totalel,1);
offset = 0;
for itk = 1:prod(npkk)
    for itx = 1:prod(npxx)
        Mat = Vcell{itk,itx};
        [MatH,MatW] = size(Mat);
        [X,Y] = meshgrid(offsetH(itk,itx)+(1:MatH),kkidx(kkIA(itk):kkIA(itk+1)-1));
        X = X';
        Y = Y';
        idx = offset+(1:MatH*MatW);
        XT(idx) = X(:);
        YT(idx) = Y(:);
        ST(idx) = Mat(:);
        if(~isempty(idx))
            offset = idx(end);
        end
    end
end
Factor.V = sparse(XT,YT,ST,totalH,totalW);

end
