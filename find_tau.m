function [best_cluster,best_tau,best_error,best_X_new] = find_tau( X, tau_range, method, n, sample_num)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

errorset = zeros(1,length(tau_range));
clusterset = zeros(n*sample_num,length(tau_range));
newdataset = cell(length(tau_range,1));
q = 1;
C = []; X_new = [];
for i = 1:length(tau_range)
    tau = tau_range(i);
    switch method
        case 'uncorrupt'
            [C,~] = lrsc( X,n,tau,q,0,0 );
        case 'noise convex'
            [C,X_new] = lrsc( X,n,tau,q,1,1 );
        case 'noise nonconvex'
            [C,X_new] = lrsc( X,n,tau,q,1,0 );
        case 'corrupt convex'
            [C,X_new] = lrsc( X,n,tau,q,2,1 );
        case 'corrupt nonconvex'
            [C,X_new] = lrsc( X,n,tau,q,2,0 );
    end
    errorset(i) = clustering_error(C,n,sample_num);
    clusterset(:,i) = C;
    newdataset{i}=X_new;
end
[best_error, error_min_index] = min(errorset);
best_tau = tau_range (error_min_index);
best_cluster = clusterset(:,error_min_index);
best_X_new = newdataset{error_min_index};

end

