function [mu,L,E] = rpca_admm(X,tau,beta,method)
%rpca_admm
% INPUT X: DxN data matrix, tau: parameter of the augmented Lagrangian,
% method: 'L1' for gross errors or 'L2' for outliers
% OUTPUT L: Low rank completion of the matrix X, E: matrix of errors

%% Define needed variables
% lambda = 1/sqrt(length(X));
N = size(X,2);
D = size(X,1);
E = zeros(size(X));
L2 = zeros(size(X));
L = zeros(size(X));
% beta = 10; 
lambda = tau;
pre_val = L + E; 
measure = 1; % initial value for convergence comparison
threshold = 10^-10;
iteration = 0;
while measure > threshold && iteration < 1000 % stop at threshold or max num. of it
    [U, S, V] = svd(X - E + L2/beta,0);
    L = U*( sign(S).*max(abs(S) - 1/beta,0) )*V';
    Z = X - L + L2/beta;
    switch method
        case 'L1'
            E = sign(Z).*max(abs(Z)-lambda/beta,0);
        case 'L21'
            z_norm = sqrt(sum(Z.*Z,1));
            zz = diag(z_norm);
            E = Z*(sign(zz).*max(abs(zz)-lambda/beta,0))/zz;
            
        otherwise
            disp('Unknown Type');
    end
    L2 = L2 + beta*(X - L - E);
    
    cur_val = L + E;
    measure = norm(pre_val-cur_val); % some measure of difference between iterations
    iteration = iteration + 1;
    pre_val = cur_val;
end
mu = mean(L,2);
end

