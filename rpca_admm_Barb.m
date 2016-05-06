function [mu,L,E] = rpca_admm(X,tau,method)
%rpca_admm
% INPUT X: DxN data matrix, tau: parameter of the augmented Lagrangian,
% method: 'L1' for gross errors or 'L2' for outliers
% OUTPUT L: Low rank completion of the matrix X, E: matrix of errors

%% Define needed variables
lambda = 1/sqrt(length(X));
N = size(X,2);
D = size(X,1);
E = zeros(size(X));
L2 = zeros(size(X));
L = zeros(size(X));

measure = 1; % initial value for convergence comparison
threshold = 10^-10;
iteration = 0;
while measure > threshold && iteration < 1000 % stop at threshold or max num. of it
    temp = L; % hold old C value
    [L, E, L2] = calc(X,E,L2,tau,lambda,method);
    % Calculate relevant convergence values
    measure = norm(temp - L)/norm(temp); % some measure of difference between iterations
    iteration = iteration + 1;
end
mu = mean(L,2);
end

