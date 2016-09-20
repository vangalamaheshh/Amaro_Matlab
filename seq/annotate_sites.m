function X = annotate_sites(X)

require_fields(X,{'chr','pos','tum_allele1','tum_allele2','ref_allele'});
chr_temp = X.chr;
                 
nx = slength(X); 
                  
for i=1:nx
  nb = setdiff([X.tum_allele1{i} X.tum_allele2{i}],X.ref_allele{i});
  if length(nb)==1
    X.newbase{i} = nb;
  elseif length(nb)==0
    X.newbase{i} = X.ref_allele{i};
  elseif length(nb)==2
    X.newbase{i} = X.tum_allele1{i};
    fprintf('Warning: Record %d is ref -> nonref1 + nonref2.  nonref1 has been arbitrarily chosen.\n');
  end
end

X = classify_muts(X);

X = rmfield(X,{'newbase','pos','ridx'});
X.chr = chr_temp;

