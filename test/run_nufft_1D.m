function run_nufft_1D(N, NG, tol, fid)

addpath('../src/');

kbox = [-N/2,N/2]';
k = -N/2:N/2-1;
kk = k(:);

xbox = [0,1]';
if(~exist(sprintf('xx_%d_nufft_1D.bin', N), 'file'))
    fprintf('Generate non-uniform distribution of x from file\n');
    xx = rand(N,1)*(N-1)/N;
    binstr = sprintf('xx_%d_nufft_1D.bin', N);
    fidxx = fopen(binstr,'w');
    string = {'DblNumMat'};
    serialize(fidxx, xx, string);
    fclose(fidxx);
else
    fprintf('Read non-uniform distribution of x from file\n');
    binstr = sprintf('xx_%d_nufft_1D.bin', N);
    fidxx = fopen(binstr,'r');
    string = {'DblNumMat'};
    xx = deserialize(fidxx, string);
end

fun = @(x,k)funFT(x,k);

f = randn(N,1) + 1i*randn(N,1);

tic;
Factor = fastBF(fun,xx,xbox,kk,kbox,NG,tol);
FactorT = toc;

tic;
yy = apply_fbf(Factor,f);
ApplyT = toc;
RunT = FactorT + ApplyT;

NC = 256;
tic;
relerr = fbf_check(N,fun,f,xx,kk,yy,NC);
Td = toc;
Td = Td*N/NC;

fprintf(fid, '------------------------------------------\n');
fprintf(fid, 'N                 : %4d\n', N);
fprintf(fid, 'Chebyshev pts     : %4d\n', NG);
fprintf(fid, 'Tolerance         : %.3e\n', tol);
fprintf(fid, 'Relative Error_2  : %.3e\n', relerr);
fprintf(fid, 'Direct Time       : %.3e s\n', Td);
fprintf(fid, 'Running Time      : %.3e mins\n', RunT/60);
fprintf(fid, 'Factorization Time: %.3e mins\n', FactorT/60);
fprintf(fid, 'Applying Time     : %.3e s\n', ApplyT);
fprintf(fid, '------------------------------------------\n\n');

end
