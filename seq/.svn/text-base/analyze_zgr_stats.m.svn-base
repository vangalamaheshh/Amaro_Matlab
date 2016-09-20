function analyze_zgr_stats(infile,outfile)
% analyze_zgr_stats(infile,outfile)
%
% creates three tables:
%
% good / bad / total                   by IGR/intron/UTR/exon
% cons / noncons / total  (all)        by ""
% cons / noncons / total  (good only)  by ""

Z = load_struct(infile);
Z = require_fields_with_convert(Z,{'terr','callablebp','nmuts','name'},{'','','','categ_name'});
Z = make_numeric(Z,{'terr','callablebp','nmuts'});

fprintf('Please ignore the following "Reconciling" error messages:\n');

R=[]; r=1;
divider2 = []; l='====================='; divider2.categ = {l}; divider2.zone = {l};
divider1 = []; l='------------'; divider1.categ = {l}; divider1.zone = {l};
z = {'IGR','intron','UTR','exon','IGR|intron|UTR|exon'};
zlabel = {'IGR','intron','UTR/pro','exon','total'};
q = {'good','bad','good|bad'}; qlabel = {'good','bad','total'};
for i=1:length(q), idx1 = grep(q{i},Z.name,1);
  for j=1:length(z), idx2 = grep(z{j},Z.name,1); idx = intersect(idx1,idx2);
    R.categ{r,1} = qlabel{i}; R.zone{r,1} = zlabel{j};
    R.Nterr(r,1) = sum(Z.terr(idx)); R.Ncov(r,1) = sum(Z.callablebp(idx));
    R.covrate(r,1) = R.Ncov(r,1)/R.Nterr(r,1); R.n_muts(r,1) = sum(Z.nmuts(idx));
    R.mutrate(r,1) = 1e6*R.n_muts(r,1)/R.Ncov(r,1);
    r=r+1;
  end, R = concat_structs_keep_all_fields({R,divider1}); r=r+1;
end
R = concat_structs_keep_all_fields({R,divider2}); r=r+1;
q = {'_cons','noncons','cons'}; qlabel = {'cons','noncons','total'};
for i=1:length(q), idx1 = grep(q{i},Z.name,1);
  for j=1:length(z), idx2 = grep(z{j},Z.name,1); idx = intersect(idx1,idx2);
    R.categ{r,1} = qlabel{i};R.zone{r,1} = zlabel{j};
    R.Nterr(r,1) = sum(Z.terr(idx));    R.Ncov(r,1) = sum(Z.callablebp(idx));
    R.covrate(r,1) = R.Ncov(r,1)/R.Nterr(r,1);    R.n_muts(r,1) = sum(Z.nmuts(idx));
    R.mutrate(r,1) = 1e6*R.n_muts(r,1)/R.Ncov(r,1);
    r=r+1;
  end, R = concat_structs_keep_all_fields({R,divider1}); r=r+1;
end
R = concat_structs_keep_all_fields({R,divider2}); r=r+1;
q = {'_cons','noncons','cons'}; qlabel = {'good cons','good noncons','good total'};
idx3 = grep('good',Z.name,1);
for i=1:length(q), idx1 = grep(q{i},Z.name,1);
  for j=1:length(z), idx2 = grep(z{j},Z.name,1); idx = intersect(idx3,intersect(idx1,idx2));
    R.categ{r,1} = qlabel{i};  R.zone{r,1} = zlabel{j};
    R.Nterr(r,1) = sum(Z.terr(idx)); R.Ncov(r,1) = sum(Z.callablebp(idx));
    R.covrate(r,1) = R.Ncov(r,1)/R.Nterr(r,1); R.n_muts(r,1) = sum(Z.nmuts(idx));
    R.mutrate(r,1) = 1e6*R.n_muts(r,1)/R.Ncov(r,1);
    r=r+1;
  end, R = concat_structs_keep_all_fields({R,divider1}); r=r+1;
end
save_struct(R,outfile);

