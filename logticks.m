function c = logticks(mn, mx)

mn(mn < 1) = 1;
minb = floor(log10(mn)); maxb = floor(log10(mx));

c = cell(maxb - minb + 1, 1);
for i = minb:maxb,
  c{i + 1} = (1:9)*(10^i);
end
c = [c{:}];

c(c < mn | c > mx) = [];
