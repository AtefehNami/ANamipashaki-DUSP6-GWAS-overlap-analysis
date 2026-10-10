%% Load DEG data for:
% 1. ADHD vs control (AN_002 vs HDFn; AN_002vsHDFn.csv)
% 2. sibling vs control (AN_001 vs HDFn; AN_001vsHDFn.csv)
function gwas_pvals = DEGanalyses_HETHOMO_GWASoverlap(which_contrast, thr_logFC, thr_FDR, num_perm)

if nargin <2
    thr_logFC = 2;
    thr_FDR = 0.05;
    num_perm = 100000;
end

DEGdata = import_DEGfile(sprintf('%s.csv', which_contrast));
DEGdata.zp = sign(DEGdata.AN_002).*-log10(DEGdata.FDR);
DEGdata.logP = -log10(DEGdata.FDR);

disorders = {'ADHD3'; 'ASD'; 'BIP3'; 'MDD4'; 'SCZ3'};
gwas_pvals = table('RowNames', disorders);


for j=1:length(disorders)

    MAGMA_data = import_MAGMAfile(sprintf('pgc%s_MAGMAdefault_genes_entrezID.txt', disorders{j}));
    fprintf('%s: BONF number: %d\n', disorders{j}, size(MAGMA_data,1));

    % evaluate overlap with disorder GWAS genes same as for drugs;
    ind_sig_genes = MAGMA_data.P<(0.05/size(MAGMA_data,1));
    sig_genes = MAGMA_data(ind_sig_genes,:);
    writetable(sig_genes, sprintf('output/MAGMAgwas_genes_%s.txt', disorders{j}), 'Delimiter', '\t')

    [~, ind_DEGgwas] = intersect(DEGdata.ensembl_gene_id, sig_genes.GENE);
    DEGgwas = DEGdata(ind_DEGgwas,1:7);
    nGWAS = size(DEGgwas,1); % number of genes to select
    fprintf('%s: number of significant GWAS genes: %d\n', disorders{j}, nGWAS);

    % number of DEG fdr<0.05 is
    data_deg_fc_fdr = DEGdata(DEGdata.FDR<thr_FDR & abs(DEGdata.AN_002)>thr_logFC,:);
    data_deg_fdr = DEGdata(DEGdata.FDR<thr_FDR,:);
    fprintf('%s: number DEGs that meet FDR<0.05 and fc>2: %d\n', disorders{j}, size(data_deg_fc_fdr,1));
    fprintf('%s: number DEGs that meet FDR<0.05: %d\n', disorders{j}, size(data_deg_fdr,1));

    % the list of genes mentioned as GWAS hits with FDR<0.05 and abs FC >2
    data_gwas_fc_fdr = DEGgwas(DEGgwas.FDR<thr_FDR & abs(DEGgwas.AN_002)>thr_logFC,:);
    data_gwas_fdr = DEGgwas(DEGgwas.FDR<thr_FDR,:);

    fprintf('%s: number of significant GWAS genes that meet FDR<0.05 and fc>2: %d\n', disorders{j}, size(data_gwas_fc_fdr,1));
    fprintf('%s: number of significant GWAS genes that meet FDR<0.05: %d\n', disorders{j}, size(data_gwas_fdr,1));

    % save output
    writetable(data_gwas_fc_fdr, sprintf('output/%s_MAGMAgwas_overlap_%d_fc_%d_fdr_%s.txt',which_contrast, thr_logFC, thr_FDR, disorders{j}), 'Delimiter', '\t')
    writetable(data_gwas_fdr, sprintf('output/%s_MAGMAgwas_overlap_%d_fdr_%s.txt',which_contrast, thr_FDR, disorders{j}), 'Delimiter', '\t')


    % get p-value using permutation testing
    gwas_num_fdr_fc = zeros(num_perm,1);
    gwas_num_fdr = zeros(num_perm,1);
    % get p-value using permutation
    % from the total list of genes, select the number of GWAS-implicated genes and count, how many fit <thr_FDR >thr_logFC
    % count how many times this number is higher than in real data

    for k=1:num_perm

        % from the total list of genes, select the number of gwas genes
        rand_inds = randsample(size(DEGdata,1), length(ind_DEGgwas));
        rand_gwas = DEGdata(rand_inds,1:7);

        % count, how many fit <thr_FDR >thr_logFC
        rand_gwas_fdr_fc = rand_gwas(rand_gwas.FDR<thr_FDR & abs(rand_gwas.AN_002)>thr_logFC,:);
        rand_gwas_fdr = rand_gwas(rand_gwas.FDR<thr_FDR,:);

        % save the number
        gwas_num_fdr_fc(k) = size(rand_gwas_fdr_fc,1);
        gwas_num_fdr(k) = size(rand_gwas_fdr,1);

    end

    %fprintf('%s: perm %d\n', disorders{j}, k);



    % count how many times this number is higher than in real data
    gwas_p_fdr_fc = sum(gwas_num_fdr_fc>size(data_gwas_fc_fdr,1))/num_perm;
    gwas_p_fdr = sum(gwas_num_fdr>size(data_gwas_fdr,1))/num_perm;

    fprintf('%s: how likely it is to find more matches FC and FDR? p=%d\n', disorders{j}, gwas_p_fdr_fc);
    fprintf('%s: how likely it is to find more matches FDR? p=%d\n', disorders{j}, gwas_p_fdr);

    % save values
    gwas_pvals.p_fdr_fc(j) = gwas_p_fdr_fc;
    gwas_pvals.p_fdr(j) = gwas_p_fdr;

    writetable(gwas_pvals, sprintf('output/%s_MAGMAgwas_overlap_pvals_%d_fc_%d_fdr.txt', which_contrast, thr_logFC, thr_FDR), 'Delimiter', '\t', 'WriteRowNames', true)

end

end
