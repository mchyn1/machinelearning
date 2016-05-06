function [ clust ] = lrsc( X,n,tau,q,type,convex )
%Outputs indicator vector for clusters
[U,S,V] = svd(X);
L = diag(S);

switch type
    case 0
        %Uncorrupted Data
        L = diag(L(L~=0));
        C = V*L*V';
    case 1
        if convex == 1
            %Data with noise: convex
            ind = find(L<=tau.^(-1/2)); %shrinkage thresholding
            L(ind) = 0;
            L(L~=0) = 1- (tau.^-1)./L.^2;
            L = diag(L);
            C = V*L*V';
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
        end
    case 2
        if convex == 1
            %Data with corrupted entries: convex
            %NEED ADMM EWWW
            
        else
            %Data with uncorrupted entries: non-convex
            %use pcp/admm for robust pca to get A
            
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

