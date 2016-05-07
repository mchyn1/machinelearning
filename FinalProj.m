%%600.6 Advanced Topics in Machine Learning
%Final Project
%%Michelle Chyn and Barbara Kim
clear all;
close all;

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
data = cell(8,1);
tau_range = [1*10^-7 1*10^-6 1*10^-5 0.0001 0.001 0.01 0.1 1 10 100 1000];

%% Adding uniform noise to (20:20:80)% of pixels
% add a loop to see if using Gaussian Noise makes a difference?

percent = .2:.2:.8;
NoisySet = cell(1,length(percent));
for i = 1:length(percent)
    X_noise=X;
    ind = randperm(a*b);
    noise = unifrnd(-60,60,[round(a*b*percent(i)),1]);
    X_noise(ind(1:round(a*b*percent(i))),:)=bsxfun(@plus,X_noise(ind(1:round(a*b*percent(i))),:),noise);
    NoisySet{i} = X_noise;
end

%% Corrupting (20:20:80)% of pixels
% 1) by removing them
% 2) by replacing with outliers of some type
for i = 1:length(percent)
    X_corrupt=X;
    ind = randperm(a*b);
    noise = unifrnd(-256,256,[round(a*b*percent(i)),1]);
    X_corrupt(ind(1:round(a*b*percent(i))),:)=bsxfun(@plus,X_corrupt(ind(1:round(a*b*percent(i))),:),noise);
    data{i+4} = X_corrupt;
end


%% LRSC with noise
data = NoisySet{2};
%convex
[C, tau, error, X_new] = find_tau(data,tau_range,'noise convex',num,c);
tau
error
figure;
scatter(1:n*c,C);
for i = 1:c*n
    figure;
    x = reshape(X(:,i),a,b);
    subplot(3,c,i);imshow(x,'DisplayRange',[]);
end

%non-convex
[C, tau, error, X_new] = find_tau(data,tau_range,'noise nonconvex',num,c);
tau
error
figure;
scatter(1:n*c,C);

%% LRSC with Uncorrupted Entries
[C, tau, error, X_new] = find_tau(data,tau_range,'uncorrupt',num,c);
tau
error
figure;
scatter(1:n*c,C);

% %% LRSC with corruptions
% %convex
% [C, tau, error] = find_tau(data,tau_range,'corrupt convex',num,c);
% tau
% error
% figure;
% scatter(1:n*c,C);
% 
% %non-convex
% [C, tau, error] = find_tau(data,tau_range,'corrupt nonconvex',num,c);
% tau
% error
% figure;
% scatter(1:n*c,C);

