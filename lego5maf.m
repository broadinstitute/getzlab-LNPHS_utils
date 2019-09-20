function [z c] = lego5maf(X,P)
% lego5maf(X,P)

pentamer_coverage_file = '/cga/tcga-gsc/home/lawrence/db/hg19/pentamer/breakdowns.mat';

if exist('P','var') && ischar(P)
  fprintf('Using "%s" as the coverage model.\n',P);
  tmp=P;
  P=[];
  P.coverage = tmp;
end

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'display','counts');
P = impose_default_value(P,'coverage','exome');

if grepmi('^(flat|counts|n)$',{P.coverage})||strcmp('counts',P.display), P.display='counts'; N = 1e6*ones(64,1); fprintf('Showing raw counts (no normalization to territory or coverage)\n');
elseif grepmi('exome|rates',{P.coverage}), load(trimer_coverage_file,'K'); N = K.Nx(1:64); fprintf('Normalizing by exome territory\n');
elseif grepmi('genome|wgs',{P.coverage}), load(trimer_coverage_file,'K'); N = K.Ng(1:64); fprintf('Normalizing by genome territory\n');
else
  error('don''t know what coverage scheme to use');
end

if ischar(X), X = load_struct(X); end

if isfield(X,'classification'), X = reorder_struct(X,strcmp('SNP',X.classification)); end
if isfield(X,'Variant_Type'), X = reorder_struct(X,strcmp('SNP',X.Variant_Type)); end
if isfield(X,'is_indel'), X = make_boolean(X, 'is_indel'); X = reorder_struct(X,~X.is_indel); end

%%%%%%%%%%%

if ~isfield(X,'context1025')
  if ~isfield(X,'pos') && isfield(X,'start'), X.pos=X.start; end
  if ~isfield(X,'pos') && isfield(X,'st'), X.pos=X.st; end
  X.context1025 = get_context1025(X.chr,X.pos);
end

if ~isfield(X,'newbase_idx')
  if ~isfield(X,'newbase')
    X.newbase = find_newbase(X);
  end
  X.newbase_idx = listmap(X.newbase,{'A','C','G','T'});
end

X = make_numeric(X,{'context1025','newbase_idx'});

midx = find(X.context1025>=1 & X.context1025<=1024 & X.newbase_idx>=1 & X.newbase_idx<=4);
try
  n = hist2d_fast(X.context1025(midx),X.newbase_idx(midx),1,1024,1,4);
catch
  n = hist2d(X.context1025(midx),X.newbase_idx(midx),1:1024,1:4);
end

%%%%%%%%%%%

load(pentamer_coverage_file,'K');   
if grepmi('exome',{P.coverage}), N = K.Nx(1:1024);
elseif grepmi('genome|wgs',{P.coverage}), N = K.Ng(1:1024);
elseif grepmi('^(flat|counts|n)$',{P.coverage}), N = 1e6*ones(1024,1);
else error('Unknown P.coverage %s',{P.coverage});
end

if isfield(X,'patient'), npat = length(unique(X.patient));
elseif isfield(X,'pat_idx'), npat = length(unique(X.pat_idx));
elseif isfield(X,'patient_idx'), npat = length(unique(X.patient_idx));
elseif isfield(X,'Tumor_Sample_Barcode'), npat = length(unique(X.Tumor_Sample_Barcode));
else npat=1; end

N=N*npat;

Nn = [N n];
Nn = collapse_Nn_1024_by_strand(Nn);
K = reorder_struct(K,1:512);

% what to use
datatype=nan;
if strcmp(P.display,'counts'), datatype=1;
elseif strcmp(P.display,'rates'), datatype=2;
else error('Unknown P.display %s',P.display); end

% format to LEGO5 grid (8x4)x(12x4)
z = nan(8*4,12*4);

base      ='ACGT';
flip(base)='TGCA';
skew(base)='CATG';
tsit(base)='GTAC';
types = {flip skew tsit};
order='TCAG';

for fromi=1:2
  from = base((3-fromi));
  by = fromi;
  for typei=1:3
    type = types{typei};
    bx = typei;
    to = type(from);
    toi = listmap(to,{'A','C','G','T'});
    for lefti=1:4, for righti=1:4
      left=order(lefti); right=order(righti);
      cx = 4*(bx-1)+righti;
      cy = 4*(by-1)+lefti;
      for flanklefti=1:4, for flankrighti=1:4
        flankleft=order(flanklefti); flankright=order(flankrighti);
        dx = 4*(cx-1)+flankrighti;
        dy = 4*(cy-1)+flanklefti;
        name = [from ' in ' flankleft left '_' right flankright];
        ki = find(strcmp(name,K.name)); if length(ki)~=1, error('?'); end
        val = Nn(ki,1+toi);
        if datatype==2, val = val / Nn(ki,1); end   % for rate, divide counts by N
        z(dy,dx) = val;
      end,end
    end,end  
  end
end

lego5(z,P);

if nargout==0
  clear z c
end
