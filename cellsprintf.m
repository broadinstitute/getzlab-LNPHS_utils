function s = cellsprintf(fmt,varargin)
% x = cellsprintf(fmt,varargin)
%
% Mike Lawrence 2010-04-22

nf = nargin-1;

x = cell(nf,1);
type = nan(nf,1);
for i=1:nf
  x{i} = varargin{i};
  if i==1
    len=length(x{i});
  else
    if length(x{i})~=len, error('Inconsistent lengths'); end
  end
  if isnumeric(x{i}) || islogical(x{i})
    type(i) = 1;
  elseif iscell(x{i}) && ischar(x{i}{1})
    type(i) = 2;
  else
    error('Unknown type');
  end
end

s = cell(len,1);
a = cell(nf,1);
for i=1:len
  for j=1:nf
    if type(j)==1, a{j} = x{j}(i);
    elseif type(j)==2, a{j} = x{j}{i};
    end
  end
  s{i} = sprintf(fmt,a{:});
end

