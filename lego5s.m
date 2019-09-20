function lego5s(varargin)
% lego5s(input1,input2,input3)
% lego5s({input1,input2,input3},titles)
% lego5s({input1,input2,input3},titles,P)


if nargin==1 && iscell(varargin{1})
  z = varargin{1};
  titles = [];
  P=[];
elseif nargin==2 && iscell(varargin{1}) && iscell(varargin{2})
  z = varargin{1};
  titles = varargin{2};
  P=[];
elseif nargin==3 && iscell(varargin{1}) && iscell(varargin{2}) && (ischar(varargin{3})&&strcmp('genome',varargin{3}))
  z = varargin{1};
  titles = varargin{2};
  P = [];
  P.coverage = 'genome';
elseif nargin==3 && iscell(varargin{1}) && iscell(varargin{2}) && (ischar(varargin{3})&&strcmp('exome',varargin{3}))
  z = varargin{1};
  titles = varargin{2};
  P = [];
  P.coverage = 'exome';
elseif nargin==3 && iscell(varargin{1}) && iscell(varargin{2}) && (ischar(varargin{3})&&strcmp('counts',varargin{3}))
  z = varargin{1};
  titles = varargin{2};
  P = [];
  P.display = 'counts';
elseif nargin==3 && iscell(varargin{1}) && iscell(varargin{2}) && isstruct(varargin{3})
  z = varargin{1};
  titles = varargin{2};
  P = varargin{3};
else
  z = varargin;
  titles = [];
  P=[];
end

P = impose_default_value(P,'nrows',[]);
P = impose_default_value(P,'fontsize',8);

n = length(z);
if ~isempty(P.nrows)
  nrows=P.nrows;
  ncols = ceil(n/nrows);
else
  if n<=3
    nrows=1;
    ncols=n;
  else
    s = ceil(sqrt(n));
    nrows=s;
    ncols=s;
    while((nrows-1)*ncols>=n) nrows=nrows-1; end
  end
end

fi=1;
for y=1:nrows, for x=1:ncols
    if fi<=n
      subplot('position',[(x-1)*(1/ncols) 1-(y*(1/nrows)) (1/ncols) (1/nrows)]);
      lego5(z{fi},P);
      if ~isempty(titles)
        zl=zlim;
        text(1,1,zl(2)*1.6,titles{fi},'fontsize',P.fontsize,'interpreter','none');
      end
    end
    fi=fi+1;
end,end,set(gcf,'color',[1 1 1]);

