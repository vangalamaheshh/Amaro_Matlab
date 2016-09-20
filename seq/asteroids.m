function asteroids(M,gene_to_analyze)

if ischar(gene_to_analyze)
  gi = find(strcmp(gene_to_analyze,M.gene.name));
elseif isnumeric(gene_to_analyze)
  % (ok)
else
  error('unknown format for gene_to_analyze');
end

nnon_g = round(sum(M.n_nonsilent(:,end,:),3));
nsil_g = round(sum(M.n_silent(:,end,:),3));
nflank_g = round(sum(M.n_flank(:,end,:),3));
Nnon_g = round(sum(M.N_non_cov(:,end,:),3));
Nsil_g = round(sum(M.N_sil_cov(:,end,:),3));
Nflank_g = round(sum(M.N_flank_cov(:,end,:),3));
globalrate_non = sum(nnon_g)/sum(Nnon_g);
globalrate_sil = sum(nsil_g)/sum(Nsil_g);
globalrate_flank = sum(nflank_g)/sum(Nflank_g);

% cosmetic changes to covariates for improved display by pr()
M.V.name = regexprep(M.V.name,'smoothed_CCLE_expression','expr');
M.V.name = regexprep(M.V.name,'rt_extra1_max','RT');
ng = M.ng;
nv = slength(M.V);
vi = find(strcmp(M.V.name,'expr'));
if ~isempty(vi)
  M.V.val{vi} = round(M.V.val{vi});
end
vi = find(strcmp(M.V.name,'hiC'));
if ~isempty(vi)
  if ~any(abs(M.V.val{vi})>1)
    M.V.val{vi} = round(1000*M.V.val{vi});
  end
end
V = nan(ng,nv);
for vi=1:nv
  V(:,vi) = M.V.val{vi};
end


% compute posterior on F, given the nsf/Nsf of this gene

gname = M.gene.name{gi};
nnon = nnon_g(gi);
nsil = nsil_g(gi);
nflank = nflank_g(gi);
Nnon = Nnon_g(gi);
Nsil = Nsil_g(gi);
Nflank = Nflank_g(gi);

F=[];
F.f = as_column([0,0.1*1.20.^[0:30]]);
tmp = binopdf(nsil,Nsil,globalrate_sil*F.f); F.post_sil = tmp/sum(tmp);
tmp = binopdf(nflank,Nflank,globalrate_flank*F.f); F.post_flank = tmp/sum(tmp);
tmp = binopdf(nnon,Nnon,globalrate_non*F.f); F.post_non = tmp/sum(tmp);

% graphical form

symb = ' .oO@';
f=F.post_sil;tmp=(f>0.0001)+(f>0.001)+(f>0.01)+(f>0.1);fh_sil = symb(tmp+1);
f=F.post_flank;tmp=(f>0.0001)+(f>0.001)+(f>0.01)+(f>0.1);fh_flank = symb(tmp+1);
f=F.post_non;tmp=(f>0.0001)+(f>0.001)+(f>0.01)+(f>0.1);fh_non = symb(tmp+1);

fprintf('\n');
fprintf('  |%s|  %-4d  %8d   %10s nonsilent    %8d  %4d %5d\n',fh_non,nnon,Nnon,gname,V(gi,1),V(gi,2),V(gi,3));
fprintf('\n');
fprintf('  |%s|  %-4d  %8d   %10s silent\n',fh_sil,nsil,Nsil,gname);
fprintf('  |%s|  %-4d  %8d   %10s flank\n',fh_flank,nflank,Nflank,gname);
fprintf('\n');

for vi=1:nv
  if isnan(V(gi,vi)) | isinf(V(gi,vi)), continue; end   % (missing value)
  dist = V(gi,vi) - V(:,vi);
  [tmp ord] = sort(abs(dist));
  numshown=0;
  for ni=1:1000
    ngi = ord(1+ni);
  
    nei_gname = M.gene.name{ngi};
    nei_nnon = nnon_g(ngi);
    nei_nsil = nsil_g(ngi);
    nei_nflank = nflank_g(ngi);
    nei_Nnon = Nnon_g(ngi);
    nei_Nsil = Nsil_g(ngi);
    nei_Nflank = Nflank_g(ngi);

    if nei_Nsil>=2e5
      tmp = binopdf(nei_nsil,nei_Nsil,globalrate_sil*F.f); F.post_sil = tmp/sum(tmp);
      symb = ' .oO@';
      f=F.post_sil;tmp=(f>0.0001)+(f>0.001)+(f>0.01)+(f>0.1);fh_sil = symb(tmp+1);
      fprintf('  |%s|  %-4d  %8d   %10s silent       %8d  %4d %5d\n',fh_sil,nei_nsil,nei_Nsil,nei_gname,V(ngi,1),V(ngi,2),V(ngi,3));
      numshown=numshown+1;
      if numshown==20, break; end
    end

    if nei_Nflank>=2e5
      tmp = binopdf(nei_nflank,nei_Nflank,globalrate_flank*F.f); F.post_flank = tmp/sum(tmp);
      symb = ' .oO@';
      f=F.post_flank;tmp=(f>0.0001)+(f>0.001)+(f>0.01)+(f>0.1);fh_flank = symb(tmp+1);
      fprintf('  |%s|  %-4d  %8d   %10s flank        %8d  %4d %5d\n',fh_flank,nei_nflank,nei_Nflank,nei_gname,V(ngi,1),V(ngi,2),V(ngi,3));
      numshown=numshown+1;
      if numshown==20, break; end
    end

  end % next neighbor

  fprintf('\n');

end % next dimension



