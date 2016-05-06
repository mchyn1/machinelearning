function [ clust ] = nspectclust( W,k )
%   INPUT: k - # of clusters, W - NxN affinity matrix
%   OUTPUT: clust - Nx1 clustering vector
d = sum(W, 2);
D = diag(d);
L = D-W;
M = diag(d.^(-1/2));
[U,S] = eig(M*L*M);
[~,i] = sort(diag(S),'ascend');
U = U(:,i(1:k));
U = normc(U);
clust = kmeans(U, k);

end

