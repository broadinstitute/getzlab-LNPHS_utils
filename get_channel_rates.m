function [ch96r ch1536r ch96f ch1536f c65 c1025] = get_channel_rates(varargin)
% get_channel_rates(<mutation struct M>, <parameters P>)
%
% returns [<96 channel rates>, <1536 channel rates>, <96 channel member>, <1536 channel member>, 
%          <context65>, <context1025>]
%
% 2014-12-05//Julian Hess

if nargin < 1, error('Not enough arguments.'); end
M = varargin{1};
if nargin == 2, P = varargin{2}; end

%parameters
%if ~exist('P'), P = []; end;
%P = impose_default_value(P, 'use_pctaf', 1);

%process fields
demand_fields(M, {'chr', 'pos'});
if ~isfield(M, 'pat') && ~isfield(M, 'patient'),
  disp('Assuming single patient MAF.');
  M.pat = repmat({'null'}, slength(M), 1);
end
if isfield(M, 'patient'), M = rename_field(M, 'patient', 'pat'); end
Mlen = slength(M);
if ~isfield(M, 'class') && ~isfield(M, 'classification'),
  disp('WARNING: assuming MAF is 100% SNPs.')
  ord = 1:slength(M);
%  error('Missing classification field.');
else
  if isfield(M, 'class'),
    [M, ord] = reorder_struct(M, strcmp(M.class, 'SNP') & cellfun('length', M.newbase) == 1);
  else
    [M, ord] = reorder_struct(M, strcmp(M.classification, 'SNP') & cellfun('length', M.newbase) == 1);
  end
end
M = make_numeric(M, {'chr', 'pos', 'context65', 'context1025'});

b('ACGTN') = 0:4;

[M, si] = sort_struct(M, {'chr', 'pos'});
[~, cui] = unique(M.chr);

if isfield(M, 'newbase'),
  M.nnewbase = cellfun(@(x) b(x), M.newbase);
elseif isfield(M, 'newbase_idx')
  disp('Assuming [ACGTN]->[0-4] newbase indexing');
  M = rename_field(M, 'newbase_idx', 'nnewbase');
end

d = NaN(slength(M), 1);
c1025n = 0; c65n = 0;
if ~isfield(M, 'context1025'), M.context1025 = d; c1025n = 1; end %3 4 2 5 1
if ~isfield(M, 'context65'), M.context65 = d; c65n = 1; end %2 1 3
M.ch1536 = d; M.ch96 = d;

for i = [cui [cui(2:end) - 1; slength(M)] M.chr(cui)]',
  j = i(1); k = i(2); l = i(3);

  if c1025n || c65n, 
    gr = upper(genome_region(l, M.pos(j:k) - 2, M.pos(j:k) + 2, 'hg19'));
    if size(gr, 1) == 1, gr = {gr}; end
  end

  %calculate context if necessary (consistent with lego5/lego3)
  if c1025n,
    M.context1025(j:k) = cellfun(@(x) sum([4^4 4^3 4^2 4 1].*b([x(3) x(4) x(2) x(5) x(1)])), gr) + 1;
  end
  if c65n,
    M.context65(j:k) = cellfun(@(x) sum([4^2 4 1].*b([x(3) x(2) x(4)])), gr) + 1;
  end

  %channels
  c1025 = double(dec2base(M.context1025(j:k) - 1, 4)) - 48;
  c65 = double(dec2base(M.context65(j:k) - 1, 4)) - 48;
  c1025 = [zeros(size(c1025, 1), 5 - size(c1025, 2)) c1025]; %pad right if necessary
  c65 = [zeros(size(c65, 1), 3 - size(c65, 2)) c65];

  nb = M.nnewbase(j:k);
  for u = as_row(find(c65(:, 1) > 1)), %take revcomp if necessary
    c65(u, [1 3 2]) = 3 - c65(u, :);
    c1025(u, [1 3 2 5 4]) = 3 - c1025(u, :);
    nb(u) = 3 - nb(u);
  end
  x = [nb c65(:, 1)];
  nbc = sum(bsxfun(@times, x(:, [2 1]), [4 1]), 2);
  nbc(nbc > 4) = nbc(nbc > 4) - 1; nbc = nbc - 1; %six mutation types, 0-5
   
  M.ch1536(j:k) = sum(bsxfun(@times, [nbc c1025(:, 2:5)], [4^4 4^3 4^2 4 1]), 2) + 1;
  M.ch96(j:k) = sum(bsxfun(@times, [nbc c65(:, 2:3)], [4^2 4 1]), 2) + 1;
end

%coverage: TODO: customize coverage profiles via P struct
pcovfile = 'ref/pentamer/breakdowns.mat';
load(pcovfile, 'K');
Kp = K;
tcovfile = 'ref/context65/breakdowns.mat';
load(tcovfile, 'K');
Kt = K;

%expand coverage to 1536/96
t1025 = dec2base(0:1023, 4) - 48;
rc1025 = sum(bsxfun(@times, 3 - t1025(:, [1 3 2 5 4]), [4^4 4^3 4^2 4 1]), 2) + 1;
T1536 = [repmat(Kp.Nx(1:256) + Kp.Nx(rc1025(1:256)), 3, 1); ...
         repmat(Kp.Nx(257:512) + Kp.Nx(rc1025(257:512)), 3, 1)];

t65 = dec2base(0:63, 4) - 48;
rc65 = sum(bsxfun(@times, 3 - t65(:, [1 3 2]), [4^2 4 1]), 2) + 1; 
T96 = [repmat(Kt.Nx(1:16) + Kt.Nx(rc65(1:16)), 3, 1); ...
       repmat(Kt.Nx(17:32) + Kt.Nx(rc65(17:32)), 3, 1)];

%names
%TODO: put this into lookup table
bc = {'A->C' 'A->G' 'A->T' 'C->A' 'C->G' 'C->T'};
B = {'A' 'C' 'G' 'T'};

tmp = dec2base((1:96) - 1, 4) - 48;
name96 = strcat(B(tmp(:, end - 1) + 1), '(', bc(sum(bsxfun(@times, tmp(:, 1:end - 2), [4 1]), 2) + 1), ')', B(tmp(:, end) + 1))';

tmp = dec2base((1:1536) - 1, 4) - 48;
name1536 = strcat(B(tmp(:, end) + 1), B(tmp(:, end - 2) + 1), '(', ...
	   bc(sum(bsxfun(@times, tmp(:, 1:end - 4), [4 1]), 2) + 1), ')', ...
	   B(tmp(:, end - 3) + 1), B(tmp(:, end - 1) + 1))';

%output channels/contexts
d = NaN(Mlen, 1);
ch96f = d; ch1536f = d;
c65 = d; c1025 = d;

%inverse sort
sii(si) = 1:slength(M);
ch96f(ord) = M.ch96(sii); ch1536f(ord) = M.ch1536(sii);
c65(ord) = M.context65(sii); c1025(ord) = M.context1025(sii);

%output rates
ch96r = [];
ch96r.rate = length(unique(M.pat))*histc(M.ch96, 1:96)./T96;
ch96r.name = as_column(name96);
ch96r.num = as_column(1:96);

ch1536r = [];
ch1536r.rate = length(unique(M.pat))*histc(M.ch1536, 1:1536)./T1536;
ch1536r.name = as_column(name1536);
ch1536r.num = as_column(1:1536); 

end
