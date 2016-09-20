function analyze_mutsig_spectra(M,outstem)

% get total mutation rates

M.patient.N_tot = sum(squeeze(M.N_cov(:,M.TOT,:)))';
M.patient.nsil_tot = sum(squeeze(M.n_silent(:,M.TOT,:)))';
M.patient.nnon_tot = sum(squeeze(M.n_nonsilent(:,M.TOT,:)))';
M.patient.rate_sil = M.patient.nsil_tot ./ M.patient.N_tot;
M.patient.rate_non = M.patient.nnon_tot ./ M.patient.N_tot;
M.patient.rate_tot = M.patient.rate_sil + M.patient.rate_non;

% get rid of hypermutated samples

avg = mean(M.patient.rate_tot);
stdev = std(M.patient.rate_tot); 
z = (M.patient.rate_tot-avg)/stdev;
thresh = 2;
remove = find(z>thresh);
M.mut.remove = ismember(M.mut.patient,remove);
M.mut.use = ismember((1:slength(M.mut))',M.use);

% get C->A, C->G, C->T, T->A, T->C, T->G rates

class = {'C->A', 'C->G', 'C->T', 'T->A', 'T->C', 'T->G'};
m = reorder_struct(M.mut,strcmp('SNP',M.mut.classification) & M.mut.use & ~M.mut.remove);
bases = {'A';'C';'G';'T'};
m.from = listmap(m.ref_allele,bases);
m.to = listmap(m.newbase,bases);
map = [0 6 5 4;1 0 2 3;3 2 0 1; 4 5 6 0];
for i=1:slength(m), m.class{i,1} = class{map(m.from(i),m.to(i))}; end
xcount(m.patient_name,m.class)

% run automatic category discovery

N = sum(M.cov.orig_cov)';
n = hist2d_fast(str2double(m.context),m.to,1,65,1,4);
Nn = collapse_Nn_65_to_32([N n]);
P=[];
P.max_k = 5;
P.mutcategs_report_filename = [outstem '.mutcateg_discovery.txt'];
Ks = find_mut_categs(Nn,P);

keyboard

