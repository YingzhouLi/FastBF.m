function y = apply_mfbf(Factors, x)
    if( size(x,2)==0 )
        y = x;
        return;
    end

    y = zeros(size(x));

    for i=1:size(Factors,1)-1
        y = y + apply_fbf(Factors{i,1},x(Factors{i,2},:));
    end
    y = y + Factors{end,1}*x(Factors{end,2},:);
end
