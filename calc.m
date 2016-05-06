function [L, E, L2] = calc(X, E, L2, tau, lambda, method)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
switch method
    case 'L1'
        L =  (X - E + L2/tau);
        [U,S,V] = svd(L);
        eps = 1/tau;
        S(S > -eps & S < eps) = 0;
        S(S > eps) = S(S > eps) - eps;
        S(S < -eps) = S(S < -eps) + eps;
        L = U*S*V';
        E = X - L + L2/tau;
        eps = lambda/tau;
        E(E >= -eps & E <= eps) = 0;
        E(E > eps) = E(E > eps) - eps;
        E(E < -eps) = E(E < -eps) + eps;
        L2 = L2 + tau*(X - L - E);
    case 'L2'
end
end

