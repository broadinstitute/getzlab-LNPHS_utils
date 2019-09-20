function pon = get_pon(chr,pos,ponfile)

% PoN returned as 8 columns:
% 1 (-) not covered (<8 total reads)
% 2 (+) covered (>=8 total reads), >99.9% ref
% 3 (1) covered, at least 1 nonref, at least 0.1% nonref
% 4 (2) covered, at least 2 nonref, at least 0.3% nonref
% 5 (3) covered, at least 3 nonref, at least 1% nonref
% 6 (4) covered, at least 3 nonref, at least 3% nonref
% 7 (5) covered, at least 3 nonref, at least 20% nonref
% 8 (6) covered, at least 10 nonref, at least 20% nonref

if ~isnumeric(chr), chr = convert_chr(chr); end
if ~isnumeric(pos), pos = str2double(pos); end
if length(chr)==1 & length(pos)>1, chr = repmat(chr,size(pos)); end
if length(chr) ~= length(pos), error('chr and pos should be same length'); end
n = length(chr);

chrlen = get_chrlen('hg19');
sz = get_filesize(ponfile);
if sz ~= 16*sum(chrlen)
  error('ponfile is not correct length for hg19');
end

coord = nan(n,1);   % 0-based
offset = 0;
for c=1:24
  idx = find(chr==c & pos>=1 & pos<=chrlen(c));
  coord(idx) = pos(idx)-1+offset;
  offset = offset + chrlen(c);
end

if ~exist(ponfile,'file'), error(['File not found: ' ponfile]); end
f = fopen(ponfile,'rb');
bof = -1;   % beginning of file

pon = nan(n,8);

[tmp ord] = sort(coord);
tt = tic;
for j=1:n, i = ord(j);
  if ~mod(j,1e5), fprintf('%d/%d ',j,n); toc(tt); end
  if isnan(coord(i)), continue; end
  fseek(f,16*coord(i),bof);
  pon(i,:) = fread(f,8,'uint16');
end
 
fclose(f);
