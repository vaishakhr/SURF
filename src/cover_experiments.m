%experiments on real data over many values of n

%the dataset is of size 58509
%code works for upto n=2^15
data = importdata('covt1.dat');
data = randsample(data,numel(data));

n = ceil(numel(data)/2);
n_ex = ceil(n^.1);

% sorting and scaling samples
samp = sort(data);
loc1 = samp(n_ex);
loc2 = samp(n+n_ex);
scale = loc2-loc1;
samp = (samp(n_ex+1:n+n_ex-1)-loc1)/scale;

gauss = fitdist(transpose(samp),'Kernel','Kernel','normal');
kern = fitdist(transpose(samp),'Kernel','Kernel','epanechnikov');

for i=1:6
    n=2^(8+i)
    n_ex = ceil(n^.1);
    samp_train = data(1:n+2*n_ex-1);
    samp_train = sort(samp_train);
    loc1 = samp_train(n_ex);
    loc2 = samp_train(n+n_ex);
    scale = loc2-loc1;
    samp_train = (samp_train(n_ex+1:n+n_ex-1)-loc1)/scale;

    % tuning parameter for SURF fixed at .25 for all our experiments
    alp = 0.25;
    deg = 1;

    [I, koi] = surf(samp_train, alp, deg);
    pres = 4*n+1;
    Y = linspace(0,1,pres+1);
    estim = @(x) regpoly(x,koi(nnz(I<=x),:));
    % the different estimators
    estim_disc = arrayfun(estim,Y);
    gauss_disc = pdf(gauss,Y);
    kern_disc = pdf(kern,Y);
    % plot(Y, estim_disc);
    % hold on;
    % plot(Y, gauss_disc);
    % plot(Y, kern_disc);
    % hold off;
    
    % precision is of the order of 1/n >> 1/sqrt(n)
    l1_diff_gauss = sum(abs(gauss_disc-estim_disc))/pres
    l1_diff_epanech = sum(abs(kern_disc-estim_disc))/pres
end