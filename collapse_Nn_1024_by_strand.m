function out = collapse_Nn_1024_by_strand(in)
% rows:  1024 categories as in /cga/tcga-gsc/home/lawrence/db/hg19/pentamer
% columns:  [N ->A ->C ->G ->T]
% Mike Lawrence 2014-11-14

if size(in,1)~=1024, error('input must have 1024 rows'); end
if size(in,2)~=5, error('input must have 5 columns (N A C G T)'); end

out = in(1:512,:,:,:,:);
compbase('ACGT') = 'TGCA';
X = generate_categ_context1025_names();
for i=1:512
  oldname = X.name{i};
  % X in XX_XX
  % 1234567890
  newname = [compbase(oldname(1)) ' in ' compbase(oldname(10)) compbase(oldname(9)) '_' compbase(oldname(7)) compbase(oldname(6))];
  j = find(strcmp(newname,X.name));
  out(i,1,:,:,:) = in(i,1,:,:,:) + in(j,1,:,:,:);
  out(i,2,:,:,:) = in(i,2,:,:,:) + in(j,5,:,:,:);
  out(i,3,:,:,:) = in(i,3,:,:,:) + in(j,4,:,:,:);
  out(i,4,:,:,:) = in(i,4,:,:,:) + in(j,3,:,:,:);
  out(i,5,:,:,:) = in(i,5,:,:,:) + in(j,2,:,:,:);
end
