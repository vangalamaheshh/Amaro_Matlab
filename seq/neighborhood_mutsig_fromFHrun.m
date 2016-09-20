function G = neighborhood_mutsig_fromFHrun(stem,P)

pat = load_struct([stem '.patients.counts_and_rates.txt']);
mut = load_struct([stem '.final_analysis_set.maf']);
load('/xchip/cga1/lawrence/mut/analysis/20110909_pancan/data.with_indels.v4.MVG.v2.mat','M','V');
G=[]; G.rank=(1:slength(M.gene))'; G.name = M.gene.name;
mut.is_sil = strcmp('Silent',mut.type) | strcmp('Synonymous',mut.type);
mut.gidx = listmap(mut.gene_name,G.name);
G.nsil = histc(mut.gidx(mut.is_sil),1:slength(G));
G.nnon = histc(mut.gidx(~mut.is_sil),1:slength(G));
patient_non_repcov = M.ttype.gene_non_repcov(1,:,end)'/M.ttype.npat(1);
patient_sil_repcov = M.ttype.gene_sil_repcov(1,:,end)'/M.ttype.npat(1);
npat = slength(pat);
G.Nnon = round(patient_non_repcov*npat);
G.Nsil = round(patient_sil_repcov*npat);

G = neighborhood_mutsig(G,V,P);
