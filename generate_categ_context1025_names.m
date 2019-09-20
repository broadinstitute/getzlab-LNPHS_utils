function categ_list = generate_categ_context1025_names()

categ_list=[];
categ_list.num = (1:1024)';
b = 'ACGT';
bi = nan(1024,5);
for i=1:1024
  t=i-1; for j=1:5, bi(i,j) = mod(t,4)+1; t=floor(t/4); end
  categ_list.name{i,1} = [b(bi(i,5)) ' in ' b(bi(i,[1 3])) '_' b(bi(i,[4 2]))];
end
categ_list.num(end+1) = 1025;
categ_list.name{end+1} = 'any N';
