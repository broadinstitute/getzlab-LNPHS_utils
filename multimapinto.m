function X = multimapinto(X,Y,xkeys,ykeys,yfields,yfields_rename)
% X = multimapinto(X,Y,xkeys[,ykeys[,yfields[,yfields_rename]]])
% Mike Lawrence 2014

if nargin<3, error('needs at least 3 args'); end

if ~exist('ykeys','var'), ykeys = xkeys; end
if ~exist('yfields','var'), yfields = fieldnames(Y); end
if exist('yfields_rename','var')
  if length(yfields)~=length(yfields_rename), error('length(yfields)~=length(yfields_rename)'); end
else
  yfields_rename = yfields;
end

if ~isstruct(X) || ~isstruct(Y), error('X and Y should be structs'); end

if ischar(xkeys), xkeys = {xkeys}; end
if ischar(ykeys), ykeys = {ykeys}; end
if ~iscellstr(xkeys) || ~iscellstr(ykeys)
  error('xkeys and ykeys should be cells-of-strings');
end

if ~all(isfield(X,xkeys)), error('not all xkeys in X'); end
if ~all(isfield(Y,ykeys)), error('not all ykeys in Y'); end

if ischar(yfields), yfields = {yfields}; end
if ~iscell(yfields)
  error('yfields should be a string or cell array of strings');
end

idx = multimap(X,Y,xkeys,ykeys);
ii = find(~isnan(idx));
if isempty(ii), fprintf('Warning: multimapinto matched no entries\n'); end

for i=1:length(yfields), fld=yfields{i};
  if ismember(fld,ykeys), continue; end  % don't overwrite the keys
  yy = getfield(Y,fld);
  if isfield(X,yfields_rename{i})
    xx = getfield(X,yfields_rename{i});
    xx(ii,:,:,:,:,:,:,:,:,:) = yy(idx(ii),:,:,:,:,:,:,:,:,:);  % don't put overwrite previous values with NaN where missing
  else
    xx = nansub(yy,idx);
  end
  X = setfield(X,yfields_rename{i},xx);
end

