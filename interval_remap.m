function c1 = interval_remap(c0, s0, e0, s1, e1) 

% Given original interval I0 = [s0 e0] and new interval I1 = [s1 e1],
% maps coordinates c0 in I0 to c1 in I1.
%
% 2016-08-22//Julian Hess

c1 = (c0 - s0)/(e0 - s0)*(e1 - s1) + s1;
