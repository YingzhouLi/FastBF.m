close all;
clear all;
%clc;

addpath('../src/');
data_path = './data/';
log_path = './log/';

if(~exist(data_path, 'dir'))
    mkdir(data_path);
end
if(~exist(log_path, 'dir'))
    mkdir(log_path);
end



% Set up parameters
i = 5;
N = 2^i;
tol = 1e-4;
NG = 2;  % number of Chebyshev pts

kbox = [-N/2,N/2-1;-N/2,N/2-1]';
k1 = rand(N,1)*(N-1)-N/2;
k2 = rand(N,1)*(N-1)-N/2;
[k1,k2] = ndgrid(k1,k2);
kk = [k1(:) k2(:)];

xbox = [0,(N-1)/N;0,(N-1)/N]';
x1 = rand(N,1)*(N-1)/N;
x2 = rand(N,1)*(N-1)/N;
[x1,x2] = ndgrid(x1,x2);
xx = [x1(:) x2(:)];

fun = @(x,k)funF(x,k);


if(1)
    f = randn(N,N) + sqrt(-1)*randn(N,N);
    binstr = sprintf('f_%d_2D.bin', N);
    fid = fopen(binstr,'w');
    string = {'CpxNumMat'};
    serialize(fid, f, string);
else
    binstr = sprintf('f_%d_2D.bin', N);
    fid = fopen(binstr,'r');
    string = {'CpxNumMat'};
    f = deserialize(fid, string);
end
f = reshape(f,N^2,1);

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

disp(['------------------------------------------']);
disp(['N                 : ' num2str(N)]);
disp(['Chebyshev pts     : ' num2str(NG)]);
disp(['Tolerance         : ' num2str(tol)]);
disp(['Relative Error_2  : ' num2str(relerr)]);
disp(['Direct Time       : ' num2str(Td) ' s']);
disp(['Running Time      : ' num2str(RunT/60) ' mins']);
disp(['Factorization Time: ' num2str(FactorT/60) ' mins']);
disp(['Applying Time     : ' num2str(ApplyT) ' s']);
disp(['------------------------------------------']);

save([data_path 'Factor_nufft_' num2str(N) '_' num2str(NG) '_2D.mat'],'Factor','-v7.3');
fid = fopen([log_path 'Factor_nufft_' num2str(N) '_' num2str(NG) '_2D.log'],'w+');
fclose(fid);
