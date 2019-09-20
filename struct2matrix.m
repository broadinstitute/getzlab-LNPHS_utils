function M = struct2matrix(S,f)
% struct2matrix(S,f)
%
% given a structure S and list of fieldnames f,
% returns a matrix M
%   in which each column in M is the (numeric) field from S
%
% if no list of fields is provided, tries to convert whole struct.
%
% Mike Lawrence 2009-03-03

if ~exist('f','var'), f = fieldnames(S); end

nf = length(f);
ns = slength(S);
M = nan(ns,nf);
for i=1:nf
  x = getfield(S,f{i});
  if ~isnumeric(x) & ~islogical(x), error('%s is neither numeric nor logical field',f{i}); end
  if islogical(x), x = double(x); end
  M(:,i) = x;
end
