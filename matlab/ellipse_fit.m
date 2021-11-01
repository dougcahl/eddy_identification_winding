function [d,l] = ellipse_fit(x,y)
M = cov(x,y);
[d,l] = eig(M); % columns are eig vecs
l = diag(l);
end
