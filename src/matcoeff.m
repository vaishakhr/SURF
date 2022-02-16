% a is coeff vector here, size of a = deg+1
% multiplication by n just for scaling
function p = matcoeff(n1,n2,n,deg)
p = zeros(1,deg+1);
for i=1:deg+1
    p(i) = n*(n2^i-n1^i)/i;
end