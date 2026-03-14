% returns coefficients that fit empvec of dimension deg+1 on interval I
function r = coeffint(samp, I, n, deg)
opnd1 = [0,1,0,0,0,0; 0,0.5,1,0,0,0 ; 0,0.2599,0.7401,1,0,0; 0,.1548,0.5,.8452,1,0; 0,0.1015,0.348,0.652,0.8985,1];
empvec = zeros(1,deg+1);
% nodes normalised for given interval
lenI = I(2)-I(1);
opnd = opnd1(deg+1,:)*lenI+I(1);
sampinI = sum(samp > I(1) & samp < I(2));
% normalization constant so that in the end, mass of 
% each interval is (#samp inside + 1)/n
normali = lenI*(sampinI+1)/(lenI*sampinI+1);
for i=1:deg+1
    % empvec is smoothed to give a hist poly for low sample intervals
%     empvec(i) = (opnd1(deg+1,i+1)-opnd1(deg+1,i) ...
%         + sum(samp > opnd(i) & samp < opnd(i+1)));
    % empvec below is extra smoothed for smaller intervals
    empvec(i) = ((opnd1(deg+1,i+1)-opnd1(deg+1,i))/lenI ...
        + sum(samp > opnd(i) & samp < opnd(i+1)))*normali*(deg+1);
end
% shifting opnd to suit present interval
matco = zeros(deg+1,deg+1);
for j = 1:deg+1
    % deg+1 is just scaling for better numerical performance,
    % rescaled back in coeff routine
    matco(j,:) = (deg+1)*matcoeff(opnd(j),opnd(j+1), n, deg);
end
% same as r = inv(matco)*transpose(empvec)
r = transpose(matco\transpose(empvec));
end