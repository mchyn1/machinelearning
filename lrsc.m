function [ clust,X_new ] = lrsc( X,n,tau,q,type,convex )
%Outputs indicator vector for clusters
%   INPUT: X - DxN data matrix, n - # of subspaces, tau, q
%   OUTPUT: clust - Nx1 data matrix,
[U,S,V] = svd(X);
L = diag(S); %vector of singular values

switch type
    case 0
        %Uncorrupted Data
        L = diag(L(L~=0));
        C = V*L*V';
        X_new = X;
    case 1
        if convex == 1
            %Data with noise: convex
            ind = L > tau.^(-1/2); %shrinkage thresholding
            V1 = V(:,ind);
            L1 = diag(L(ind).^-2);
            C = V1*(eye(size(L1,1)) - (1/tau)*L1)*V1.';
            X_new = X*C;
        else
            %Data with noise: non-convex
            ind = find(abs(L)<=sqrt(2./tau));
            L(ind) = 0;
            l = diag(L);
            L = zeros(size(S));
            L(1:size(l,1),1:size(l,2)) = l;
            A = U*L*V';
            [~,Sa,Va] = svd(A);
            La = diag(Sa);
            La = diag(La(La~=0));
            C = Va*La*Va';
            X_new = A;
        end
    case 2
        if convex == 1
            %Data with corrupted entries: convex
            B = size(X,1)*size(X,2)/(4*norm(X,1));
            [C,~] = admm(X,B,convex);
            X_new = X*C;
        else
            %Data with uncorrupted entries: non-convex
            %use pcp/admm for robust pca to get A
            B = size(X,1)*size(X,2)/(4*norm(X,1));
            [A,~] = admm(X,B,convex);
            X_new = A;
            
            %using theorem 8.6 on A to ge tC
            [~,Sa,Va] = svd(A);
            La = diag(Sa);
            La = diag(La(La~=0));
            C = Va*La*Va';
        end
end

W = abs(C).^q;
clust = nspectclust(W,n);

end

