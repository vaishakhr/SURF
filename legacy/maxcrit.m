% recursive function to compute, verified using I with 8*[1/8] and [1/2,1/8,1/8,1/8,1/8]
function r = maxcrit(samp, I, cumuprobI, maincoeff, n, deg, alp)
coeffI = coeffint(samp, [I(1),I(end)], n, deg);
sumprob = cumuprobI(numel(cumuprobI));
% criteria evaluated on the given I
differI = differ(I, maincoeff, coeffI)...
    -alp*sqrt(5*(deg+1)*sumprob*log(n)/n);
if numel(I) == 2
    % if there is just one interval do 
    r = differI;
else
    % separating into intervals if there is more than one interval in I
    I1 = I(cumuprobI<=.5*sumprob);
    I2 = I(cumuprobI>=.5*sumprob);
    cumuprobI1 = cumuprobI(cumuprobI<=.5*sumprob);
    cumuprobI2 = cumuprobI(cumuprobI>=.5*sumprob)-.5*sumprob;
    samp1 = samp(samp >= I1(1) & samp <= I1(2));
    samp2 = samp(samp >= I2(1) & samp <= I2(2));
    r = max(differI,...
        maxcrit(samp1, I1, cumuprobI1, maincoeff, n, deg, alp)+...
        maxcrit(samp2, I2, cumuprobI2, maincoeff, n, deg, alp));
end