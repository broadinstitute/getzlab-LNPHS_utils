function tier = get_tiers(gene,tiersfile)

if ~exist('tiersfile','var'), tiersfile = '/cga/tcga-gsc/home/lawrence/mut/20140218_pancan/gene_tiers.v9.mat'; end

load(tiersfile,'G');

tier = mapacross(standardize_mutsig2cv_genenames(gene),standardize_mutsig2cv_genenames(G.gene),G.tier);
