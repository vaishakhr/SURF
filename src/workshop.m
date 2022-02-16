% gaussian
% degree
deg = 1;
% tuning parameter for SURF
alp = 0.25;
% mixture probability
prob = .65;
% number of samples
n = 2^11;
% some extra samples tossed out, penalty added in l_1 error
n_ex = ceil(n^0.1);

l1_diff = 0;
for i = 1:10
    m = n+2*n_ex;
    f = rand(1,m) <= prob;
    f = floor(floor(2*f)/2);
    % gaussian mixture
    gauss_m = -.45;
    gauss_sd = .15;
    gauss_m_1 = .3;
    gauss_sd_1 = .2;
    samp = f.*normrnd(gauss_m, gauss_sd, [1,m])+...
    (1-f).*normrnd(gauss_m_1, gauss_sd_1, [1,m]);
    % remove 2n_ex samples out and rescale
    samp = sort(samp);
    loc1 = samp(n_ex);
    loc2 = samp(n+n_ex);
    scale = loc2-loc1;
    samp = (samp(n_ex+1:n+n_ex-1)-loc1)/scale;

    [I,koi] = surf(samp, alp, deg);

%     disp("Number of pieces in the partition is")
%     numel(I)

    % now onto plotting
    % precision to calculate l_1 difference
    pres = 16*n+1;
    Y = linspace(0,1,pres+1);
    % the distribution and estimator respectively
    gaussdensity = @(x) prob*exp(-(x-gauss_m)^2/(2*gauss_sd^2))/(2*pi*gauss_sd^2)^0.5...
         +(1-prob)*exp(-(x-gauss_m_1)^2/(2*gauss_sd_1^2))/(2*pi*gauss_sd_1^2)^0.5;
    estim = @(x) regpoly(x,koi(nnz(I<=x),:));
    density_disc1 = scale*arrayfun(gaussdensity,Y*scale+loc1);
    estim_disc = arrayfun(estim,Y);

%     % plotting
%     hold off;
%     plot(Y,density_disc1);
%     hold on;
%     plot(Y,estim_disc);

    % l1_diff = integral(@(x) absdiff(x),0,1);
    % discrete l1 difference, 32 just to take care of blow up at 0,
    % overall discretized error order is 1/n << 1/\sqrt{n}, the learning error
    l1_diff = l1_diff+sum(abs(density_disc1(32:end)-estim_disc(32:end)))/pres+2*n_ex/n;
end
l1_diff/10