log_path = './log/';

if(~exist(log_path, 'dir'))
    mkdir(log_path);
end

NGlist = [4 6 8];
tollist = [1e-5 1e-7 1e-10];

func_list = {'fun0'};
for func_i = 1:length(func_list)
    func_name = func_list{func_i};
    for N = 2.^(5:2:9)
        fid = fopen([log_path 'Factor_nufft_2D_' num2str(N) '.log'],'a+');
        for iter = 1:length(NGlist)
            NG = NGlist(iter);
            tol = tollist(iter);
            run_nufft_2D(N, NG, tol, fid);
            fprintf('Func %s, N %4d, NG %2d finished.\n',func_list{func_i},N,NG);
        end
        fclose(fid);
    end
end
