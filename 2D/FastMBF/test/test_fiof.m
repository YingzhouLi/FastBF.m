close all;
clear all;
clc;

addpath('../src/');
data_path = './data/';
log_path = './log/';

if(~exist(data_path, 'dir'))
    mkdir(data_path);
end
if(~exist(log_path, 'dir'))
    mkdir(log_path);
end


for i=6:2:6
%% Set up parameters
N = 2^i; %powers of 2 from 64 to 65536;
tol = 1e-4;
NG = 5;  % number of Chebyshev pts

kbox = [-N/2,N/2-1;-N/2,N/2-1];
k = -N/2:N/2-1;
[k1,k2] = ndgrid(k);
kk = [k1(:) k2(:)];

xbox = [0,(N-1)/N;0,(N-1)/N];
x = (0:N-1)/N;
[x1,x2] = ndgrid(x);
xx = [x1(:) x2(:)];

func_name = 'funF';
switch func_name
case 'funF'
    fun = @(x,k)funF(x,k);
case 'fun0'
    fun = @(x,k)fun0(x,k);
case 'fun1'
    fun = @(x,k)fun1(x,k);
case 'fun2'
    fun = @(x,k)fun2(x,k);
end


%% Begin test
if(1)
    f = randn(N,N) + sqrt(-1)*randn(N,N);
    binstr = sprintf('f_%d.bin', N);
    fid = fopen(binstr,'w');
    string = {'CpxNumMat'};
    serialize(fid, f, string);
end
if(0)
    binstr = sprintf('f_%d.bin', N);
    fid = fopen(binstr,'r');
    string = {'CpxNumMat'};
    f = deserialize(fid, string);
end
f = reshape(f,N^2,1);

tic;
Factors = fastMBF(N,fun,xx,xbox,kk,kbox,NG,tol);
FactorT = toc;

tic;
yy = apply_mfbf(Factors,f);
ApplyT = toc;
RunT = FactorT + ApplyT;

NC = 256;
tic;
relerr = fbf_check(N,fun,f,yy,NC);
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

%save([data_path 'Factor_' func_name '_' num2str(N) '_' num2str(NG) '.mat'],'Factor','-v7.3');
%fid = fopen([log_path 'Factor_' func_name '_' num2str(N) '_' num2str(NG) '.log'],'w+');
%fclose(fid);

end
