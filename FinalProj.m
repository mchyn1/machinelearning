%%600.6 Advanced Topics in Machine Learning
%Final Project
%%Michelle Chyn and Barbara Kim

%% loading images
num = 3; %num of subspaces (individuals)
a = 50; %192
b = 50; %168
c = 6; %10 %images per individual
s = [a, b, c];
X = [];
%X = zeros(a*b,c);
for n = 1:num %number of individuals/subspaces
    x = zeros(192,168,10);
    for i = 1:c %for each image from 1 individual
        x(:,:,i) = loadimage(n,i);
    end
    x = x(1:a,1:b,1:c);
    x = reshape(x,[a*b,c]);
    X = [X x];
end

%% Adding uniform noise to (20:20:80)% of pixels
% add a loop to see if using Gaussian Noise makes a difference?
percent = .2:.2:.8;
for i = 1:length(percent)
    X_noise=X;
    ind = randperm(a*b);
    noise = unifrnd(-256,256,[round(a*b*percent(i)),1]);
    X_noise(ind(1:round(a*b*percent(i))),:)=bsxfun(@plus,X_noise(ind(1:round(a*b*percent(i))),:),noise);
end

%% Corrupting (20:20:80)% of pixels
% 1) by removing them
% 2) by replacing with outliers of some type


%% LRSC with noise

%convex
tau = 1*10^-6; %.00001;
q = 1;

[C,X_new] = lrsc( X_noise,n,tau,q,1,1 );
error = clustering_error(C,num,c)
figure;
scatter(1:n*c,C);

%non-convex
tau = 10000; %.00001;
q = 1;
runs = 10;
for r = 1:runs
    [C,X_new] = lrsc( X_noise,n,tau,q,1,0 );
    figure;
    scatter(1:n*c,C);
end


%% LRSC with corruptions
%convex

%non-convex


%% LRSC with Uncorrupted Entries
tau = 10000; %.00001;
q = 1;
runs = 10;
for r = 1:runs
    [C,~] = lrsc( X,n,tau,q,0,0 );
    figure;
    scatter(1:n*c,C);
end





