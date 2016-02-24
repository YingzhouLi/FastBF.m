log_path = './log/';

if(~exist(log_path, 'dir'))
    mkdir(log_path);
end

func_list = {'fun0'};
for func_i = 1:length(func_list)
    func_name = func_list{func_i};
    for N = 2.^(8:2:18)
        fid = fopen([log_path 'Factor_' func_name '_' num2str(N) '.log'],'a+');
        for NG = 4:2:8
            tol = 1e-4;
            run_fastbf_1D(N, func_name, NG, tol, fid);
            fprintf('Func %s, N %4d, NG %2d finished.\n',func_list{func_i},N,NG);
        end
        fclose(fid);
    end
end
