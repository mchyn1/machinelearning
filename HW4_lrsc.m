%% Problem 4, LRSC
clear all;
close all;

s = [62,48]; %size of image (original size is [192,168])
sample = 10; %10 sample images

dim = s(1)*s(2); % dimension of image
n = 3; % number of subspaces
imgset = zeros(dim,n*sample);
q = 1;
for individual_number = 1:n % individual number
    for i = 1:sample
        image = loadimage(individual_number,i);
        imgset(:,q) = reshape(image(1:s(1),1:s(2)), [dim 1]);
        q = q + 1;
    end
end

%%
tau_range = [1*10^-7 1*10^-6 1*10^-5 0.0001 0.001 0.01 0.1 1 10 100];
errorset = zeros(1,length(tau_range)); clusterset = zeros(n*sample,length(tau_range));
for p = 1:length(tau_range)
    tau = tau_range(p);
    q = 1;
    cluster = lrsc(imgset,n,tau,q,1,1);
    ord = 1:n*sample;
    used_subspace = [];
    error_mat = zeros(n,1);
    for i = 1:n
        k = 1;
        for j = 1:n
            temp = cluster(k:k+(sample-1));
            boop(j) = sum(temp == i);
            k = k + sample;
        end
        [values, max_ind] = sort(boop,'descend');
        m = 1;
        while any(used_subspace == max_ind(m))
            m = m + 1;
        end
        max_value = values(m);
        used_subspace = [used_subspace max_ind(m)];
        error_mat(i) = sample - max_value;
    end
    error = sum(error_mat)/(n*sample);
    errorset(p) = error; % contains errors for each tau
    clusterset(:,p) = cluster; % contains clusters for each tau
end
%% Find best tau and error
[error_min_value, error_min_index] = min(errorset)
final_tau = tau_range (error_min_index)
scatter(ord,clusterset(:,error_min_index));