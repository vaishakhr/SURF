%data preprocessing from salary.text
%the dataset is of size 32561
%code works for upto n=2^14
sal = importdata('salary.dat');

% degree
deg = 2;
% number of samples 
n = 2^12;
n_ex = ceil(n^0.5);
sal_train = sal(1:n+2*n_ex-1);
sal_test = sal(n:n+n/4-1);
% tuning parameter for SURF
alp = 0.25;

% removing outliers, sorting and scaling samples
samp = sort(sal_train);
loc1 = samp(n_ex);
loc2 = samp(n+n_ex);
scale = loc2-loc1;
samp = (samp(n_ex+1:n+n_ex-1)-loc1)/scale;

[I, koi] = surf(samp, alp, deg);
gauss = fitdist(transpose(samp),'Normal');
kern = fitdist(transpose(samp),'Kernel','Kernel','epanechnikov');
beta = fitdist(transpose(samp),'Beta');
logist = fitdist(transpose(samp),'Logistic');
% pres = 16*n+1;
% Y = linspace(0,1,pres+1);
% % the estimator respectively
% estim = @(x) regpoly(x,koi(nnz(I<=x),:));
% estim_disc = arrayfun(estim,Y);
% gauss_disc = pdf(gauss,Y);
% kern_disc = pdf(kern,Y);
% plot(Y, estim_disc);
% hold on;
% plot(Y, gauss_disc);
% plot(Y, kern_disc);
% hold off;

m = n/4;
m_ex = ceil(m^0.5);
samp_test =  sal(n+2*n_ex:n+2*n_ex+m+2*m_ex);
samp_test = sort(samp_test);
loc1 = samp_test(m_ex);
loc2 = samp_test(m+m_ex);
scale = loc2-loc1;
samp_test = (samp_test(m_ex+1:m+m_ex-1)-loc1)/scale;
% unfair to all estimators to consider outliers since 
% none learns in KL distance, we trim outliers
% to compare
samp_test_trim = samp_test(m_ex:m-m_ex);

likeli_SURF = arrayfun(estim,samp_test_trim);
perplex_SURF = exp(-sum(log(likeli_SURF))/numel(samp_test_trim))

likeli_gauss = pdf(gauss,samp_test(m_ex:m-m_ex));
perplex_gauss = exp(-sum(log(likeli_gauss))/numel(samp_test_trim))

likeli_kern = pdf(kern,samp_test(m_ex:m-m_ex));
perplex_kern = exp(-sum(log(likeli_kern))/numel(samp_test_trim))

likeli_beta = pdf(beta,samp_test(m_ex:m-m_ex));
perplex_beta = exp(-sum(log(likeli_beta))/numel(samp_test_trim))

likeli_logist = pdf(beta,samp_test(m_ex:m-m_ex));
perplex_logist = exp(-sum(log(likeli_logist))/numel(samp_test_trim))
