%%600.6 Advanced Topics in Machine Learning
%Final Project
%%Michelle Chyn and Barbara Kim
clear all;
close all;

%% loading images
num = 3; %num of subspaces (individuals)
a = 50; %192
b = 50; %168
c = 10; %10 %images per individual
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
data = cell(9,1);
data{1} = X;
tau_range = [1*10^-7 1*10^-6 1*10^-5 0.0001 0.001 0.01 0.1 1 10 60 100 150 1000];

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
CorruptSet = cell(1,length(percent));
for i = 1:length(percent)
    X_corrupt=X;
    ind = randperm(a*b);
    noise = unifrnd(-256,256,[round(a*b*percent(i)),1]);
    X_corrupt(ind(1:round(a*b*percent(i))),:)=bsxfun(@plus,X_corrupt(ind(1:round(a*b*percent(i))),:),noise);
    CorruptSet{i} = X_corrupt;
end


%% LRSC with noise
for i = 1:4
    data{i+1} = NoisySet{i};
    data{i+5} = CorruptSet{i};
end
method = cell(5,1);
method{1} = 'uncorrupt';
method{2} = 'noise convex';
method{3} = 'noise nonconvex';
method{4} = 'corrupt convex';
method{5} = 'corrupt nonconvex';
method_tau = zero(d,test);
method_error = zeros(d,test);
method_computeTime = zeros(d,test);
for d = 3 %1:9
    for test = 1:3 %1:5
        tic
        [C, tau, error, X_new] = find_tau(data{d},tau_range,method{test},num,c);
        method_computeTime = toc;
        method_tau(d,test) = tau;
        method_error(d,test) = error;
        figure;
        scatter(1:n*c,C);
        figure;
        for i = 1:c*n
            if test ==1
                X_new = data{d};
            end
            x = reshape(X_new(:,i),a,b);
            subplot(n,c,i);imshow(x,'DisplayRange',[]);
        end
    end
end

%{
%non-convex
[C, tau, error, X_new] = find_tau(data,tau_range,'noise nonconvex',num,c);
tau
error
figure;
scatter(1:n*c,C);
figure;
for i = 1:c*n
    x = reshape(X_new(:,i),a,b);
    subplot(n,c,i);imshow(x,'DisplayRange',[]);
end

%% LRSC with Uncorrupted Entries
% [C, tau, error, X_new] = find_tau(data,tau_range,'uncorrupt',num,c);
% tau
% error
% figure;
% scatter(1:n*c,C);

% %% LRSC with corruptions
% %convex
% [C, tau, error, X_new] = find_tau(data,tau_range,'corrupt convex',num,c);
% tau
% error
% figure;
% scatter(1:n*c,C);
%
% %non-convex
% [C, tau, error, X_new] = find_tau(data,tau_range,'corrupt nonconvex',num,c);
% tau
% error
% figure;
% scatter(1:n*c,C);

%}