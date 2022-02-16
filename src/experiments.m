% this code contains synthetic experiments that were run using SURF
% parameters alp, deg, n
% n must be a power of 2
% code written for deg <= 4
% a nice value of the parameter alp = 0.25
% code is not MATLAB optimized and written recursively,
% running time may be large depending on n, deg

% % beta mixture
% % degree
% deg = 1;
% % tuning parameter for SURF
% alp = 0.25;
% % beta mixture
% beta_alp = .8;
% beta_bet = 4;
% beta_alp_1 = 2; 
% beta_bet_1 = 2;
% % mixture probability
% prob = .4;
% % number of samples
% n = 2^9;

% l1_diff = 0;
% for i=1:10
%     f = rand(1,n-1) <= prob;
%     f = floor(floor(2*f)/2);
%     samp = f.*betarnd(beta_alp, beta_bet, [1,n-1])+...
%     (1-f).*betarnd(beta_alp_1, beta_bet_1, [1, n-1]);
%     samp = sort(samp);
%     
%     [I,koi] = surf(samp, alp, deg);
% 
%     disp("Number of pieces in the partition is")
%     numel(I)
% 
%     now onto plotting
%     precision to calculate l_1 difference
%     pres = 16*n+1;
%     Y = linspace(0,1,pres+1);
%     the distribution and estimator respectively
%     betadensity = @(x) prob*(x^(beta_alp-1))*((1-x)^(beta_bet-1))/beta(beta_alp,beta_bet)...
%          +(1-prob)*(x^(beta_alp_1-1))*((1-x)^(beta_bet_1-1))/beta(beta_alp_1,beta_bet_1);
%     estim = @(x) regpoly(x,koi(nnz(I<=x),:));
%     density_disc1 = arrayfun(betadensity,Y);
%     estim_disc = arrayfun(estim,Y);
% 
%     l1_diff = integral(@(x) absdiff(x),0,1);
%     discrete l1 difference, 32 just to take care of blow up at 0,
%     overall discretized error order is 1/n << 1/\sqrt{n}, the learning error
%     l1_diff = l1_diff+sum(abs(density_disc1(32:end)-estim_disc(32:end)))/pres;
% end
% l1_diff/10

% % plotting the last run
% hold off;
% plot(Y,density_disc1);
% hold on;
% plot(Y,estim_disc);


%-----------------%-----------------%-----------------

% gamma
% degree
deg = 1;
% tuning parameter for SURF
alp = 0.25;
% gamma mixture
% mixture probability
prob = .2;
gam_a = 4;
gam_b = .04;
gam_a_1 = 8;
gam_b_1 = 0.06;
% number of samples
n = 2^10;
% some extra samples tossed out to take care of tails, 
% penalty added in l_1 error
n_ex = ceil(n^0.1);

l1_diff = 0;
for i=1:10
    m = n+n_ex;
    f = rand(1,m) <= prob;
    f = floor(floor(2*f)/2);
    samp = f.*gamrnd(gam_a, gam_b, [1,m])+...
    (1-f).*gamrnd(gam_a_1, gam_b_1, [1,m]);
    % remove n_ex samples out on tail and rescale
    samp = sort(samp);
    scale = samp(n);
    samp = samp(1:n-1)/scale;

    [I,koi] = surf(samp, alp, deg);

%     disp("Number of pieces in the partition is")
%     numel(I)

    % precision to calculate l_1 difference
    pres = 16*n+1;
    Y = linspace(0,1,pres+1);
    % the distribution and estimator respectively
    gammadensity = @(x) prob*(x^(gam_a-1))*(exp(-x/gam_b))/(gamma(gam_a)*gam_b^gam_a)...
         +(1-prob)*(x^(gam_a_1-1))*(exp(-x/gam_b_1))/(gamma(gam_a_1)*gam_b_1^gam_a_1);
    estim = @(x) regpoly(x,koi(nnz(I<=x),:));
    density_disc1 = scale*arrayfun(gammadensity,Y*scale);
    estim_disc = arrayfun(estim,Y);

    % l1_diff = integral(@(x) absdiff(x),0,1);
    % discrete l1 difference, 32 just to take care of blow up at 0,
    % overall discretized error order is 1/n << 1/\sqrt{n}, the learning error
    l1_diff = l1_diff+sum(abs(density_disc1(32:end)-estim_disc(32:end)))/pres+n_ex/n;
end
l1_diff/10

% plotting the last run
hold off;
plot(Y,density_disc1);
hold on;
plot(Y,estim_disc);

% % gaussian
% % degree
% deg = 1;
% % tuning parameter for SURF
% alp = 0.25;
% % gaussian mixture
% gauss_m = .4;
% gauss_sd = .1;
% gauss_m_1 = .6;
% gauss_sd_1 = .2;
% % mixture probability
% prob = .3;
% % number of samples
% n = 2^15;
% % some outliers tossed out, penalty added in l_1 error
% n_ex = ceil(n^0.1);
% 
% l1_diff = 0;
% for i = 1:5
%     m = n+2*n_ex;
%     f = rand(1,m) <= prob;
%     f = floor(floor(2*f)/2);
%     samp = f.*normrnd(gauss_m, gauss_sd, [1,m])+...
%     (1-f).*normrnd(gauss_m_1, gauss_sd_1, [1,m]);
%     % remove 2n_ex samples out and rescale
%     samp = sort(samp);
%     loc1 = samp(n_ex);
%     loc2 = samp(n+n_ex);
%     scale = loc2-loc1;
%     samp = (samp(n_ex+1:n+n_ex-1)-loc1)/scale;
% 
%     [I,koi] = surf(samp, alp, deg);
% 
% %     disp("Number of pieces in the partition is")
% %     numel(I)
% 
%     % precision to calculate l_1 difference
%     pres = 16*n+1;
%     Y = linspace(0,1,pres+1);
%     % the distribution and estimator respectively
%     gaussdensity = @(x) prob*exp(-(x-gauss_m)^2/(2*gauss_sd^2))/(2*pi*gauss_sd^2)^0.5...
%          +(1-prob)*exp(-(x-gauss_m_1)^2/(2*gauss_sd_1^2))/(2*pi*gauss_sd_1^2)^0.5;
%     estim = @(x) regpoly(x,koi(nnz(I<=x),:));
%     density_disc1 = scale*arrayfun(gaussdensity,Y*scale+loc1);
%     estim_disc = arrayfun(estim,Y);

%     % l1_diff = integral(@(x) absdiff(x),0,1);
%     % discrete l1 difference, 32 just to take care of blow up at 0,
%     % overall discretized error order is 1/n << 1/\sqrt{n}, the learning error
%     l1_diff = l1_diff+sum(abs(density_disc1(32:end)-estim_disc(32:end)))/pres+2*n_ex/n;
% end
% l1_diff/5
% % plotting the last run
% hold off;
% plot(Y,density_disc1);
% hold on;
% plot(Y,estim_disc);