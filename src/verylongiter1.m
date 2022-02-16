% gamma
% tuning parameter for SURF
alp = 0.25;
% mixture probability
prob = .4;
gam_a = 3;
gam_b = .05;
gam_a_1 = 6;
gam_b_1 = 0.075;
gam_gam_a = factorial(gam_a-1);
gam_gam_a_1 = factorial(gam_a_1-1);

% degree
deg = 1
% number of samples
for j= 9:13
n = 2^j
% some outliers tossed out, penalty added in l_1 error
n_ex = ceil(n^0.1);

l1_diff = 0;
niter = 5*ceil(14/j);
for i=1:10
    m = n+n_ex;
    f = rand(1,m) <= prob;
    f = floor(floor(2*f)/2);
    % gamma mixture generation
    samp = f.*gamrnd(gam_a, gam_b, [1,m])+...
    (1-f).*gamrnd(gam_a_1, gam_b_1, [1,m]);
    % remove n_ex samples out on tail and rescale
    samp = sort(samp);
    scale = samp(n);
    samp = samp(1:n-1)/scale;

    [I,koi] = surf(samp, alp, deg);

%     disp("Number of pieces in the partition is")
%     numel(I)

    % now onto plotting
    % precision to calculate l_1 difference
    pres = 16*n+1;
    Y = linspace(0,1,pres+1);
    % the distribution and estimator respectively
    gammadensity = @(x) prob*(x^(gam_a-1))*(exp(-x/gam_b))/(gam_gam_a*gam_b^gam_a)...
         +(1-prob)*(x^(gam_a_1-1))*(exp(-x/gam_b_1))/(gam_gam_a_1*gam_b_1^gam_a_1);
    estim = @(x) regpoly(x,koi(nnz(I<=x),:));
    density_disc1 = scale*arrayfun(gammadensity,Y*scale);
    estim_disc = arrayfun(estim,Y);

%     % plotting
%     hold off;
%     plot(Y,density_disc1);
%     hold on;
%     plot(Y,estim_disc);

    % l1_diff = integral(@(x) absdiff(x),0,1);
    % discrete l1 difference, 32 just to take care of blow up at 0,
    % overall discretized error order is 1/n << 1/\sqrt{n}, the learning error
    l1_diff = l1_diff+sum(abs(density_disc1(32:end)-estim_disc(32:end)))/pres+n_ex/n;
end
l1_diff/niter
end


% gamma
% degree
deg = 2
% tuning parameter for SURF
alp = 0.25;
% mixture probability
prob = .2;
% number of samples
for j= 9:13
n = 2^j
% some outliers tossed out, penalty added in l_1 error
n_ex = ceil(n^0.1);

l1_diff = 0;
niter = 5*ceil(14/j);
for i=1:10
    m = n+n_ex;
    f = rand(1,m) <= prob;
    f = floor(floor(2*f)/2);
    % gamma mixture generation
    samp = f.*gamrnd(gam_a, gam_b, [1,m])+...
    (1-f).*gamrnd(gam_a_1, gam_b_1, [1,m]);
    % remove n_ex samples out on tail and rescale
    samp = sort(samp);
    scale = samp(n);
    samp = samp(1:n-1)/scale;

    [I,koi] = surf(samp, alp, deg);

%     disp("Number of pieces in the partition is")
%     numel(I)

    % now onto plotting
    % precision to calculate l_1 difference
    pres = 16*n+1;
    Y = linspace(0,1,pres+1);
    % the distribution and estimator respectively
    gam_gam_a = factorial(gam_a-1);
    gam_gam_a_1 = factorial(gam_a_1-1);
    gammadensity = @(x) prob*(x^(gam_a-1))*(exp(-x/gam_b))/(gam_gam_a*gam_b^gam_a)...
         +(1-prob)*(x^(gam_a_1-1))*(exp(-x/gam_b_1))/(gam_gam_a_1*gam_b_1^gam_a_1);
    estim = @(x) regpoly(x,koi(nnz(I<=x),:));
    density_disc1 = scale*arrayfun(gammadensity,Y*scale);
    estim_disc = arrayfun(estim,Y);

%     % plotting
%     hold off;
%     plot(Y,density_disc1);
%     hold on;
%     plot(Y,estim_disc);

    % l1_diff = integral(@(x) absdiff(x),0,1);
    % discrete l1 difference, 32 just to take care of blow up at 0,
    % overall discretized error order is 1/n << 1/\sqrt{n}, the learning error
    l1_diff = l1_diff+sum(abs(density_disc1(32:end)-estim_disc(32:end)))/pres+n_ex/n;
end
l1_diff/niter
end


% gamma
% degree
deg = 3
% tuning parameter for SURF
alp = 0.25;
% mixture probability
prob = .2;
% number of samples
for j= 9:13
n = 2^j
% some outliers tossed out, penalty added in l_1 error
n_ex = ceil(n^0.1);

l1_diff = 0;
niter = 5*ceil(14/j);
for i=1:10
    m = n+n_ex;
    f = rand(1,m) <= prob;
    f = floor(floor(2*f)/2);
    % gamma mixture generation
    samp = f.*gamrnd(gam_a, gam_b, [1,m])+...
    (1-f).*gamrnd(gam_a_1, gam_b_1, [1,m]);
    % remove n_ex samples out on tail and rescale
    samp = sort(samp);
    scale = samp(n);
    samp = samp(1:n-1)/scale;

    [I,koi] = surf(samp, alp, deg);

%     disp("Number of pieces in the partition is")
%     numel(I)

    % now onto plotting
    % precision to calculate l_1 difference
    pres = 16*n+1;
    Y = linspace(0,1,pres+1);
    % the distribution and estimator respectively
    gammadensity = @(x) prob*(x^(gam_a-1))*(exp(-x/gam_b))/(gam_gam_a*gam_b^gam_a)...
         +(1-prob)*(x^(gam_a_1-1))*(exp(-x/gam_b_1))/(gam_gam_a_1*gam_b_1^gam_a_1);
    estim = @(x) regpoly(x,koi(nnz(I<=x),:));
    density_disc1 = scale*arrayfun(gammadensity,Y*scale);
    estim_disc = arrayfun(estim,Y);

%     % plotting
%     hold off;
%     plot(Y,density_disc1);
%     hold on;
%     plot(Y,estim_disc);

    % l1_diff = integral(@(x) absdiff(x),0,1);
    % discrete l1 difference, 32 just to take care of blow up at 0,
    % overall discretized error order is 1/n << 1/\sqrt{n}, the learning error
    l1_diff = l1_diff+sum(abs(density_disc1(32:end)-estim_disc(32:end)))/pres+n_ex/n;
end
l1_diff/niter
end