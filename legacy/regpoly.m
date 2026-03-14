% this is just a regular poly, size of a is deg+1
function p=regpoly(x,a)
p = 0;
for i = 1:size(a,2)
    p = p+x.^(i-1)*a(i);
end
end