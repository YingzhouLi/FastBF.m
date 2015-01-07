function idx = fbf_subidx(zz,z1box,z2box,npz1,npz2)
    zidx = cell(1,1);
    zidx{1} = 1:size(zz,1);
    npzz1 = 1;
    npzz2 = 1;
    z1len = z1box(2)-z1box(1);
    z2len = z2box(2)-z2box(1);
    while(npzz1<npz1 | npzz2<npz2)
        if(npzz1<npz1)
            zidxcur = cell(npzz1*2,npzz2);
            for z1 = 1:npzz1
            for z2 = 1:npzz2
                tmpidx = zidx{z1,z2};
                
                z1s = z1box(1)+z1len/npzz1*(z1-1);
                z1e = z1box(1)+z1len/npzz1/2*(2*z1-1);
                
                zidxcur{2*z1-1,z2} = tmpidx(z1s<=zz(tmpidx,1) & zz(tmpidx,1)<z1e);
                
                z1s = z1e;
                if(z1==npzz1)
                    z1e = z1box(1)+z1len/npzz1*z1+1;
                else
                    z1e = z1box(1)+z1len/npzz1*z1;
                end
                zidxcur{2*z1,z2} = tmpidx(z1s<=zz(tmpidx,1) & zz(tmpidx,1)<z1e);
            end
            end
            zidx = zidxcur;
            npzz1 = npzz1*2;
        end
        if(npzz2<npz2)
            zidxcur = cell(npzz1,npzz2*2);
            for z1 = 1:npzz1
            for z2 = 1:npzz2
                tmpidx = zidx{z1,z2};
                
                z2s = z2box(1)+z2len/npzz2*(z2-1);
                z2e = z2box(1)+z2len/npzz2/2*(2*z2-1);
                
                zidxcur{z1,2*z2-1} = tmpidx(z2s<=zz(tmpidx,2) & zz(tmpidx,2)<z2e);
                
                z2s = z2e;
                if(z2==npzz2)
                    z2e = z2box(1)+z2len/npzz2*z2+1;
                else
                    z2e = z2box(1)+z2len/npzz2*z2;
                end
                zidxcur{z1,2*z2} = tmpidx(z2s<=zz(tmpidx,2) & zz(tmpidx,2)<z2e);
            end
            end
            zidx = zidxcur;
            npzz2 = npzz2*2;
        end
    end
    idx = zidx;
end