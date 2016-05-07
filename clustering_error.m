function [error] = clustering_error(cluster,n,sample)
%   INPUT: cluster - Nx1, n - # of subspaces, sample - # of images per
%   subspace
ord = 1:n*sample;
used_subspace = [];
error_mat = zeros(1,n);
for i = 1:n
    boop = zeros(1,n);
    k = 1;
    for j = 1:n
        temp = cluster(k:k+(sample-1));
        boop(j) = sum(temp == i); %sum how many points are clustered into subspace i
        k = k + sample;
    end
    %identify which subpsace the most points corresponding to this individual were clustered into
    [values, indices] = sort(boop,'descend');  
    m = 1;
    while any(used_subspace == indices(m)) %if subspace has already been used, select next subspace
        m = m + 1;
    end
    used_subspace = [used_subspace indices(m)];
    max_value = values(m); %how many points were clustered in this subspace?
    error_mat(i) = sample - max_value; %what is the clustering error for this subspace
end
error = sum(error_mat)/(n*sample);

end

