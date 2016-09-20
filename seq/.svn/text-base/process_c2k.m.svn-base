function [M G mutrate] = process_c2k(direc)

d = dir([direc '/TCGA-*']);
M = {}; np=0;
for i=1:length(d)
  fname = [direc '/' d(i).name '/mutation_reports5f.maf'];
  if exist(fname,'file')
    np=np+1;
    M{np} = load_struct(fname,'%s%s%f%f%s%s%s%s%s',0);
  end
end
M = concat_structs(M(:));
M = rename_field(M,{'col1','col2','col3','col4','col5','col6','col7','col8','col9'},...
   {'build','chr','start','end','ref_allele','tum_allele1','tum_allele2','tumor_barcode','normal_barcode'});
chr_temp = M.chr;
M.chr = convert_chr(M.chr);
M.pos = M.start;
M.newbase = M.tum_allele2;
nx = slength(M);
for i=1:nx
  nb = setdiff([M.tum_allele1{i} M.tum_allele2{i}],M.ref_allele{i});
  if length(nb)==1
    M.newbase{i} = nb;
  elseif length(nb)==0
    M.newbase{i} = M.ref_allele{i};
  elseif length(nb)==2
    M.newbase{i} = M.tum_allele1{i};
    fprintf('Warning: Record %d is ref -> nonref1 + nonref2.  nonref1 has been arbitrarily chosen.\n');
  end
end

M = classify_muts(M);
M.type = fillblanks(M.type,'IGR');

M = rmfield(M,{'newbase','pos','ridx'});
M.chr = chr_temp;

M2 = reorder_struct(M,ismember(M.type,{'Missense','Nonsense','Splice_site'}));

G = []; [G.name gi gj] = unique(M2.gene);
for i=1:slength(G), G.n(i,1) = sum(gj==i); end

S = load_struct('/xchip/tcga/gbm/analysis/lawrence/refseq/hg18/gene_length.txt','%s%s%f');
sidx = listmap(G.name,S.gene);
i2 = find(~isnan(sidx));
G.N = nan(slength(G),1);
G.N(i2) = S.exons_length(sidx(i2));
G = reorder_struct(G,~isnan(G.N));

G.r = G.n ./ G.N;

ntot = sum(G.n);
Ntot = sum(G.N);
mutrate = ntot / (Ntot * np);

G.p = 1-binocdf(G.n-1,G.N,mutrate);
G.q = calc_fdr_value(G.p);

G = sort_struct(G,'q');
