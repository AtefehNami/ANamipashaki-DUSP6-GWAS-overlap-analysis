% ============================================================
% Run GWAS overlap permutations: simple, length-matched, length and detectability matched; 
% ============================================================

thr_logFC = 2; 
thr_FDR = 0.05; 
num_perm = 100000; 

gwas_pvals_AN_002vsHOM = DEGanalyses_HETHOMO_GWASoverlap('AN_002vsHOM', thr_logFC, thr_FDR, num_perm); 
gwas_pvals_AN_002vsHET = DEGanalyses_HETHOMO_GWASoverlap('AN_002vsHET', thr_logFC, thr_FDR, num_perm); 
gwas_pvals_HETvsHOM = DEGanalyses_HETHOMO_GWASoverlap('HETvsHOM', thr_logFC, thr_FDR, num_perm); 
