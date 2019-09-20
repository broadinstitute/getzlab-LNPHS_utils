function m = load_matrix(infile,truncate_flag)
%
% Loads a textfile that contains only numbers.
% May have any number of rows.  No header row.
% Each row should have the same number of columns.
% Columns may be tab- or space-delimited.
% Multiple tabs/spaces are treated as single.
%
% Mike Lawrence 2009-06-24

if ~exist(infile,'file'), error('%s not found',infile); end

X = load_textfile(infile);

regexprep(X,'\s*','\t');   % collapse multiple delimiters

while true   % remove blank lines
  if isempty(X), break; end
  if X(end)~=char(10), break; end
  X = X(1:end-1);
end  

if isempty(X)
  m = [];
else
  nrows = sum(X==char(10))+1;
  d = sscanf(X,'%f');
  ncells = length(d);
  ncols = ncells / nrows;
  if ncols ~= round(ncols)
    if ~exist('truncate_flag','var') | ~truncate_flag
      error('Matrix contains %d rows and %d cells -> %f columns ?!',nrows,ncells,ncols);
    else
      fprintf('Truncating last row.\n');
      ncols = round(ncols);
      nrows = floor(ncells / ncols);    
    end
  end
  m = reshape(d(1:ncols*nrows),ncols,nrows)';
end
