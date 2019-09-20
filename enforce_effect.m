function eff = enforce_effect(m,P)
% returns eff = effect as in mutation_type_dictionary, one of:
%               ncd, syn, mis, non, spl, indel_ncd, indel_cod, indel_spl
%               or unk (unknown) for DNPs etc.

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'context_and_effect_fwb_file',    'ref/gencode/c65e29/all.fwb');
P = impose_default_value(P,'context_and_effect_categs_file', 'ref/gencode/c65e29/categs.txt');

if ~isfield(m,'pos') && isfield(m,'start'), m.pos=m.start; end
demand_fields(m,{'chr','pos','newbase'});
if ~isnumeric(m.chr), m.chr=convert_chr(m.chr); end
if ~isnumeric(m.pos), m.pos=str2doubleq_wrapper(m.pos); end

m.context = get_from_fwb(m.chr,m.pos,P.context_and_effect_fwb_file);

K = load_struct(P.context_and_effect_categs_file);
K.name = regexprep(K.name,'splice-site','spl/spl/spl');
K.name = regexprep(K.name,'noncoding','ncd/ncd/ncd');
K = parsein(K,'name','(.) in (.)_(.):(...)/(...)/(...)$',{'from','left','right','eff1','eff2','eff3'});
z=repmat({'---'},slength(K),1);K.effA=z;K.effC=z;K.effG=z;K.effT=z;
idx=find(strcmp(K.from,'A'));K.effC(idx)=K.eff1(idx);K.effG(idx)=K.eff2(idx);K.effT(idx)=K.eff3(idx);
idx=find(strcmp(K.from,'C'));K.effA(idx)=K.eff1(idx);K.effG(idx)=K.eff2(idx);K.effT(idx)=K.eff3(idx);
idx=find(strcmp(K.from,'G'));K.effA(idx)=K.eff1(idx);K.effC(idx)=K.eff2(idx);K.effT(idx)=K.eff3(idx);
idx=find(strcmp(K.from,'T'));K.effA(idx)=K.eff1(idx);K.effC(idx)=K.eff2(idx);K.effG(idx)=K.eff3(idx);

eff = repmat({'---'},slength(m),1);
idx = find(strcmp(m.newbase,'A'));eff(idx)=nansub(K.effA,m.context(idx));
idx = find(strcmp(m.newbase,'C'));eff(idx)=nansub(K.effC,m.context(idx));
idx = find(strcmp(m.newbase,'G'));eff(idx)=nansub(K.effG,m.context(idx));
idx = find(strcmp(m.newbase,'T'));eff(idx)=nansub(K.effT,m.context(idx));

% deal with non-SNPs
is_nonsnp = false(slength(m),1);
is_indel = false(slength(m),1);
if isfield(m,'classification')
  is_nonsnp(~strcmp('SNP',m.classification))=true;
  is_indel(grepm('INS|DEL',m.classification))=true;  
end
if isfield(m,'newbase')
  is_nonsnp(cellfun('length',m.newbase)~=1)=true;
  is_indel(grepm('\-',m.newbase))=true;
end
if isfield(m,'ref_allele')
  is_nonsnp(cellfun('length',m.ref_allele)~=1)=true;
  is_indel(grepm('\-',m.ref_allele))=true;
end
if isfield(m,'start') && isfield(m,'end')
  is_nonsnp(m.end~=m.start)=true;
end

is_nonsnp_and_nonindel = (is_nonsnp & ~is_indel);
idx = find(is_nonsnp_and_nonindel);
if ~isempty(idx)
  fprintf('%d mutations have unknown effect.  possible DNPs etc?\n',length(idx));
  eff(idx) = repmat({'unk'},length(idx),1);
end

% indels
snp_eff = nansub(K.effA,m.context); idx = find(strcmp('---',snp_eff)); snp_eff(idx) = nansub(K.effC,m.context(idx));
idx = find(is_indel & grepm('ncd',snp_eff)); eff(idx) = repmat({'indel_ncd'},length(idx),1);
idx = find(is_indel & grepm('syn|mis|non',snp_eff)); eff(idx) = repmat({'indel_cod'},length(idx),1);
idx = find(is_indel & grepm('spl',snp_eff)); eff(idx) = repmat({'indel_spl'},length(idx),1);














