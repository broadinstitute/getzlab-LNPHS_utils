function varargout = get_chrlen(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = load_chrlen(varargin{:});
else
  [varargout{1}] = load_chrlen(varargin{:});
end
