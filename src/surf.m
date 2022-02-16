% surf takes sorted samples strictly in [0,1] as input 
%(num samples must be power of 2-1),
% tuning parameter alpha (a nice value of alp = 0.25),
% degree d <= 4
% surf return piecewise intervals I in [0,1] and coefficients koi
% Note:
% code is written recursively and not MATLAB optimized and
% may take a long time to run depending on platform, number of samples, deg
function [I,koi] = surf(samp, alp, deg)
hold off;
warning('off');
n = numel(samp)+1;

% I is interval list, consists of start and end points, if merged, smaller list
I = [0, samp, 1];
% cumulative problist
cumuprobI = linspace(0, 1, n+1);
% coefficients initilaization
koi = zeros(n+1,deg+1);
% kois last element is dummy element
for i=1:n
    koi(i,:)=coeffint(samp, [I(i),I(i+1)], n, deg);
end
D = log2(n);

for i = 1:D
    for j = 1:2^(D-i)
        % Icand is interval candidate
        probinit = (j-1)*2^(-1*(D-i));
        probfin = j*2^(-1*(D-i));
        cumuprobIcand = cumuprobI(cumuprobI>=probinit & cumuprobI<=probfin)-probinit;
        Icand = I(cumuprobI>=probinit & cumuprobI<=probfin);
        [res,koooi] = merge(samp, Icand, cumuprobIcand, n, deg, alp);
        if res == 1
            I = I(cumuprobI<=probinit | cumuprobI>=probfin);
            koi = koi(cumuprobI<=probinit | cumuprobI>=probfin,:);
            cumuprobI = cumuprobI(cumuprobI<=probinit | cumuprobI>=probfin);
            koi(cumuprobI==probinit,:) = koooi;
        end    
    end
end
end

