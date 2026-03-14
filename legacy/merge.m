% returns decision to merge, I is considered interval, probI is probability
% partiton, samp are samples restricted to I, deg is degree, alp is tuning
% parameter
function [r, koi] = merge(samp, I, cumuprobI, n, deg, alp)
% this is the merged polynomial that is passed along
maincoeff = coeffint(samp, [I(1),I(end)], n, deg);
% merge will be default called with more than one candidate interval
% they're evaluated below
sumprob = cumuprobI(numel(cumuprobI));
I1 = I(cumuprobI<=.5*sumprob);
I2 = I(cumuprobI>=.5*sumprob);
cumuprobI1 = cumuprobI(cumuprobI<=.5*sumprob);
cumuprobI2 = cumuprobI(cumuprobI>=.5*sumprob)-.5*sumprob;
samp1 = samp(samp >= I1(1) & samp <= I1(2));
samp2 = samp(samp >= I2(1) & samp <= I2(2));
val = maxcrit(samp1, I1, cumuprobI1, maincoeff, n, deg, alp)+...
    maxcrit(samp2, I2, cumuprobI2, maincoeff, n, deg, alp);
% disp('merged')
if val <= 0
    r = 1;
else
    r = 0;
end
koi = maincoeff;
end