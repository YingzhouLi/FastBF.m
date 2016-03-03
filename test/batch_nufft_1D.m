log_path = './log/';

if(~exist(log_path, 'dir'))
    mkdir(log_path);
end

NGlist = [4 6 8];
tollist = [1e-4 1e-6 1e-8];

for N = 2.^(8:2:20)
    fid = fopen([log_path 'Factor_' func_name '_1D_' num2str(N) '.log'],'a+');
    for iter = 1:length(NGlist)
        NG = NGlist(iter);
        tol = tollist(iter);
        run_fastbf_1D(N, func_name, NG, tol, fid);
        fprintf('Nufft, N %4d, NG %2d finished.\n',N,NG);
    end
    fclose(fid);
end
