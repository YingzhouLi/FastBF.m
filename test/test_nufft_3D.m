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
i = 4;
N = 2^i;
tol = 1e-9;
NG = 3;  % number of Chebyshev pts

kbox = [-N/2,N/2-1;-N/2,N/2-1;-N/2,N/2-1]';
if(0)
    kk = rand(N^3,3)*(N-1)-N/2;
else
    k = -N/2:N/2-1;
    [k1,k2,k3] = ndgrid(k);
    kk = [k1(:) k2(:) k3(:)];
end

xbox = [0,(N-1)/N;0,(N-1)/N;0,(N-1)/N]';
if(1)
    xx = rand(N^3,3)*(N-1)/N;
else
    x = (0:N-1)/N;
    [x1,x2,x3] = ndgrid(x);
    xx = [x1(:) x2(:) x3(:)];
end

fun = @(x,k)funFT(x,k);

f = randn(N,N,N) + 1i*randn(N,N,N);
f = reshape(f,N^3,1);

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

%save([data_path 'Factor_nufft_' num2str(N) '_' num2str(NG) '_3D.mat'],'Factor','-v7.3');
fid = fopen([log_path 'Factor_nufft_' num2str(N) '_' num2str(NG) '_3D.log'],'w+');
fclose(fid);

