%importing data from sensor1.dat
%the dataset is of size 57923 after removing outliers, and comes shuffled
%code works for upto n=2^14, d<=2
sen = importdata('sensor1.dat');

% perplex_SURF = 0;
% perplex_gauss = 0;
% perplex_gamma = 0;
% perplex_beta = 0;
% perplex_kern_epa = 0;
% perplex_kern_gauss = 0;

% degree
deg = 1;
% number of samples 
n = 2^14;

num_rounds = floor(4/5*numel(sen)/n)-1;
niter = 0;
perplex_SURF = 0;
perplex_gauss = 0;
perplex_gamma = 0;
perplex_beta = 0;
perplex_kern_epa = 0;
perplex_kern_gauss = 0;
offset = 1;
for i=1:num_rounds
sen_train = sen(offset+1:offset+n+1);
sen_test = sen(offset+n+2:offset+n+2+n/4);
offset = offset+n+2+n/4;
% tuning parameter for SURF
alp = 0.25;

% sorting and scaling samples
samp = sort(sen_train);
loc1 = samp(1);
loc2 = samp(end);
scale = loc2-loc1;
samp = (samp(2:end-1)-loc1)/scale;

[I, koi] = surf(samp, alp, deg);
gauss = fitdist(transpose(samp),'Normal');
gamma = fitdist(transpose(samp),'Gamma');
beta = fitdist(transpose(samp),'Beta');
logist = fitdist(transpose(samp),'Logistic');
kern_gau = fitdist(transpose(samp),'Kernel','Kernel','normal');
kern_epa = fitdist(transpose(samp),'Kernel','Kernel','epanechnikov');

% pres = 16*n+1;
% Y = linspace(0,1,pres+1);
% % the estimator respectively
% estim = @(x) regpoly(x,koi(nnz(I<=x),:));
% 
% estim_disc = arrayfun(estim,Y);
% kern_epa_disc = pdf(kern_epa,Y);
% kern_gauss_disc = pdf(kern_gau,Y);
% plot(Y, estim_disc);
% hold on;
% plot(Y, kern_epa_disc);
% plot(Y, kern_gauss_disc);
% hold off;

samp_test = sort(sen_test);
loc1 = samp_test(1);
loc2 = samp_test(end);
scale = loc2-loc1;
samp_test = (samp_test(2:end-1)-loc1)/scale;

likeli_SURF = arrayfun(estim,samp_test);
if(any(likeli_SURF<0))
    continue;
end
perplex_SURF = perplex_SURF + exp(-sum(log(likeli_SURF))/numel(samp_test));

likeli_gauss = pdf(gauss,samp_test);
perplex_gauss = perplex_gauss+ exp(-sum(log(likeli_gauss))/numel(samp_test));

likeli_gamma = pdf(gamma,samp_test);
perplex_gamma = perplex_gamma+ exp(-sum(log(likeli_gamma))/numel(samp_test));

likeli_beta = pdf(beta,samp_test);
perplex_beta =  perplex_beta + exp(-sum(log(likeli_beta))/numel(samp_test));

likeli_kern_epa = pdf(kern_epa,samp_test);
perplex_kern_epa = perplex_kern_epa + exp(-sum(log(likeli_kern_epa))/numel(samp_test));

likeli_kern_gau = pdf(kern_gau,samp_test);
perplex_kern_gauss = perplex_kern_gauss+ exp(-sum(log(likeli_kern_gau))/numel(samp_test));

niter = niter + 1;
end
perplex_SURF/niter
% perplex_gauss/niter
% perplex_gamma/niter
% perplex_beta/niter
perplex_kern_epa/niter
perplex_kern_gauss/niter
