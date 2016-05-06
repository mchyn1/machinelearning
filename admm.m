function [ C,E ] = admm( X,B )

E = 0;
Lambda = 0; temp = Lambda;
lambda = 1./sqrt(max(size(X,1),size(X,2)));
tau = lambda*B^-1;
r = 0;
thresh = .0001;
while norm(temp-Lambda)/norm(temp)<thresh && r<1000
    r = r+1;
    temp = Lambda;
    arg1 = X - E + B^-1*Lambda;
    [U,S,V] = svd(arg1);
    S(S<=B^-1&&S>=-B^-1) = 0;
    S(S>B^-1) = S - B^-1;
    S(S<B^-1) = S + B^-1;
    C = U*S*V';
    arg2 = X - X*C + B^-1*Lambda;
    arg2(arg2<=tau&&arg2>=-tau) = 0;
    arg2(arg2>tau) = arg2 - tau;
    arg2(arg2<tau) = arg2 + tau;
    E = arg2;
    Lambda = Lambda + B.*(X-X*C-E);
end


end

