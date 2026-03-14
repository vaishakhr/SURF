% computes absolute differences between c1 and c2 in I
function r = differ(I, c1, c2)
r = integral(@(x) abs(regpoly(x,c1-c2)),I(1),I(2));
end
