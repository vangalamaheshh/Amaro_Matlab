function T = annotate_M_calls(M)

is_germline = (size(M,2)==9);

bases = {'A';'C';'G';'T';'N'};
models = {'AA';'AC';'AG';'AT';'CC';'GG';'TT';'CG';'CT';'GT'};

T = [];
T.chr = M(:,1);
T.pos = M(:,2);
T.ref_allele = bases(M(:,3));

if is_germline % results from germcall

  m = models(M(:,4));
  T.newbase = T.ref_allele;
  T.tum_allele1 = T.ref_allele;
  T.tum_allele2 = T.ref_allele;
  for i=1:slength(T)
    T.tum_allele1{i} = m{i}(1);
    T.tum_allele2{i} = m{i}(2);
    isref = (m{i}==T.ref_allele{i});
    if isref(1) & ~isref(2)
      T.newbase{i} = T.tum_allele2{i};
    elseif isref(2) & ~isref(1)
      T.newbase{i} = T.tum_allele1{i};
    elseif ~isref(1) & ~isref(2)
      T.newbase{i} = m{i};  % try both
    else % ref+ref (i.e. no mutation)
      T.newbase{i} = T.ref_allele{i};
    end
  end
  T.LOD = M(:,5);
  T.fMAPQZ = M(:,6);
  T.avgnMM = M(:,7);
  T.nUWP = M(:,8);
  T.filter = M(:,9);

else % results from somcall

  T.tum_allele1 = T.ref_allele;
  T.tum_allele2 = bases(M(:,4));
  T.newbase = T.tum_allele2;
  T.t_LOD = M(:,5);
  T.n_LOD = M(:,6);
  T.t_fMAPQZ = M(:,7);
  T.n_fMAPQZ = M(:,8);
  T.t_avgnMM = M(:,9);
  T.n_avgnMM = M(:,10);
  T.t_nUWP = M(:,11);
  T.n_nUWP = M(:,12);
  T.filter = M(:,13);
end

rs = compare_to_dbSNP(T.chr,T.pos);
z = find(rs==0);
nz = setdiff(1:length(rs),z);
T.dbSNP = repmat({'none'},slength(T),1);
for i=1:length(nz), T.dbSNP{nz(i)} = ['rs' num2str(rs(nz(i)))]; end
T = classify_muts(T);
T = rmfield(T,{'newbase','ridx'});
T.gene_longname = get_longnames(T.gene);
