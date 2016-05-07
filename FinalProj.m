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
percent = .2:.2:.8;
levels = .1:.1:.6;
data = cell(length(percent)*length(levels)+length(percent)+1,1);
data{1} = X;
tau_range = [1*10^-7 1*10^-6 1*10^-5 0.0001 0.001 0.01 0.1 1 10 60 100 150 1000];

%% Adding uniform noise to (20:20:80)% of pixels
% add a loop to see if using Gaussian Noise makes a difference?

NoisySet = cell(1,length(percent)*length(levels));
k = 0;
for j = 1:length(levels)
    for i = 1:length(percent)
        k = k + 1;
        X_noise=X;
        ind = randperm(a*b);
        noise = unifrnd(-255*levels(j),255*levels(i),[round(a*b*percent(i)),1]);
        X_noise(ind(1:round(a*b*percent(i))),:)=bsxfun(@plus,X_noise(ind(1:round(a*b*percent(i))),:),noise);
        NoisySet{k} = abs(X_noise);
    end
end

%% Corrupting (20:20:80)% of pixels
CorruptSet = cell(1,length(percent));
for i = 1:length(percent)
    X_corrupt=X;
    ind = randperm(a*b);
    noise = unifrnd(-255,255,[round(a*b*percent(i)),1]);
    X_corrupt(ind(1:round(a*b*percent(i))),:)=bsxfun(@plus,X_corrupt(ind(1:round(a*b*percent(i))),:),noise);
    CorruptSet{i} = X_corrupt;
end


%% LRSC
j = 1;
for i = 1:length(percent) + length(percent)*length(levels)
    if i <= length(NoisySet)
        data{i+1} = NoisySet{i};
    else
        data{i+1} = CorruptSet{j};
        j = j + 1;
    end
end
data_length = length(data);
method = cell(5,1);
method{1} = 'uncorrupt';
method{2} = 'noise convex';
method{3} = 'noise nonconvex';
method{4} = 'corrupt convex';
method{5} = 'corrupt nonconvex';
method_tau = zeros(data_length,5);
method_error = zeros(data_length,5);
method_computeTime = zeros(data_length,5);
for d = 1:24 %1:9
    for test = 1:3 %1:5
        tic
        [C, tau, error, X_new] = find_tau(data{d},tau_range,method{test},num,c);
        method_computeTime(d,test) = toc;
        method_tau(d,test) = tau;
        method_error(d,test) = error;
%         figure;
%         scatter(1:n*c,C);
%         figure;
%         for i = 1:c*n
%             if test ==1
%                 X_new = data{d};
%             end
%             x = reshape(X_new(:,i),a,b);
%             subplot(n,c,i);imshow(x,'DisplayRange',[]);
%         end
    end
end

%% Plot Figures
% Noisy Data
ordered_error = zeros(length(levels),length(percent),5);
ordered_tau = zeros(length(levels),length(percent),5);
ordered_computeTime = zeros(length(levels),length(percent),5);
for i = 1:length(percent)
    figure();
    hold on;
    xlabel('Noise Level');
    ylabel('Clustering Error (%)');
    title(sprintf('Noise in %.f%% of the Data',percent(i)*100));
    for j = 1:test %change by number of methods run
        ind = 1;
        for k = i:length(percent):length(levels)*length(percent)
            ordered_error(ind,i,j) = method_error(k,j);
            ordered_tau(ind,i,j) = method_tau(k,j);
            ordered_computeTime(ind,i,j) = method_computeTime(k,j);
            ind = ind + 1;
        end
        plot(levels,ordered_error(:,i,j)*100,'-o');
    end
    legend('Uncorrupt LRSC','Noisy Convex LRSC','Noisy Non-Convex LRSC');
end
