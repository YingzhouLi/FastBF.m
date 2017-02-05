function cidxs = childidx(siz,idx)

n = length(siz);
cidxs = zeros(1,2^n);
for itc = 1:2^n
    vec_child = idx2vec(2*ones(1,n),itc);
    vec = idx2vec(siz,idx);
    cidxs(itc) = vec2idx(2*siz, (vec-1)*2 + vec_child);
end

end