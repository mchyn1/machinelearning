function [C] = corrupt_convex_admm(X,tau,beta)
%rpca_admm
% INPUT X: DxN data matrix, tau: parameter of the augmented Lagrangian,

%% Define needed variables
% lambda = 1/sqrt(length(X));
N = size(X,2);
D = size(X,1);
E = zeros(size(X));
Z = zeros(N);
C = zeros(N);
Y1 = zeros(size(X)); Y2 = zeros(N);
% beta = 10; 
lambda = tau;
pre_val = C;
measure = 10; % initial value for convergence comparison
threshold = 10^-10;
iteration = 0;
while (measure > threshold || norm(C) == 0) && iteration < 1000 % stop at threshold or max num. of it
    [U, S, V] = svd(Z + Y2/beta,0);
    C = U*( sign(S).*max(abs(S) - 1/beta,0) )*V';
    Z = (eye(N) + X.'*X)\(X.'*(X - E) + C + (X.'*Y1 - Y2)/beta);
    Loo = X - X*Z + Y1/beta;
    E = sign(Loo).*max(abs(Loo)-lambda/beta,0);
    Y1 = Y1 + beta*(X - X*Z - E);
    Y2 = Y2 + beta*(Z - C);

    cur_val = C;
    measure = norm(pre_val-cur_val); % some measure of difference between iterations
    iteration = iteration + 1;
    pre_val = cur_val;
end
iteration
end

