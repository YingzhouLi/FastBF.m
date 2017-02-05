function pidx = parentidx(siz,idx)

vec = idx2vec(siz,idx);
pidx = vec2idx(siz/2, floor((vec+1)/2));

end