function g = standardize_mutsig2cv_genenames(g)

g = regexprep(g,'^GRLF1$','ARHGAP35');
g = regexprep(g,'^JUB$','AJUBA');
g = regexprep(g,'^MLL$' ,'KMT2A');
g = regexprep(g,'^MLL2$','KMT2D');
g = regexprep(g,'^MLL3$','KMT2C');
g = regexprep(g,'^MLL4$','KMT2B');
g = regexprep(g,'^MLL5$','KMT2E');
g = regexprep(g,'^BHD$','FLCN');
g = regexprep(g,'^COPEB$','KLF6');
g = regexprep(g,'^GRAF$','ARHGAP26');
g = regexprep(g,'^HRPT2$','CDC73');
g = regexprep(g,'^MADH4$','SMAD4');
g = regexprep(g,'^NBS1$','NBN');
g = regexprep(g,'^PTCH$','PTCH1');
g = regexprep(g,'^SDH5$','SDHAF2');
g = regexprep(g,'^TCF1$','HNF1A');
g = regexprep(g,'^TNFRSF6$','FAS');
g = regexprep(g,'^WTX$','FAM123B');
g = regexprep(g,'^SRSF2$','SFRS2');



