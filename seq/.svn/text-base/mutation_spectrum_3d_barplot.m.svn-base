function mutation_spectrum_3d_barplot(mut)

if ~isfield(mut,'newbase'), mut.newbase = find_newbase(mut); end
if ~isfield(mut,'context65')
  b = interpret_build(q.NCBI_Build{1});
  fwb = ['/xchip/cga1/lawrence/db/hg' num2str(b) '/context65/all.fwb'];
  q.context65 = get_from_fwb(q.Chromosome,q.Start_position,fwb);
else
  mut = make_numeric(mut,'context65');
end

n = hist2d_fast(mut.context65,listmap(mut.newbase,{'A';'C';'G';'T'}),1,65,1,4);

% load mean coverage from a typical set of well-covered samples (EXOMES)
load('/xchip/cga1/lawrence/dlbcl/analysis/20111003/mayo/maf_file_source/An_DLBCL_20110815.coverage.mat','C1');
N = round(collapse_categories_to_65(mean(C1.orig_cov)','/xchip/cga1/lawrence/db/hg19/c65e/categs.txt'));

% -- the value of npat doesn't actually really matter
% -- bars will be displayed relative to an invisible scale
% -- if you want to display absolute rates, change npat to reflect
%    how many patients are represented in the list of mutations
npat = 1;
N=N*npat;
Nn = collapse_Nn_65_to_32([N sum(n,3)]);

draw_mutation_spectrum_3d_barplot(Nn);
