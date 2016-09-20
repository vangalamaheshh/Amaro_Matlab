function [sets,sets_lst,sets_idx]=gsea_get_set(set_name,chipname,gacc)

gsea_dir='/home/radon00/gadgetz/gsea/';
fname=[chipname '_' set_name];


if exist([ gsea_dir fname '.mat'],'file')
  load([ gsea_dir fname '.mat']);
else
  if exist([ gsea_dir fname '.gmx'],'file')
    extname='.gmx';
  else
    extname='.gmt';
  end
  sets=read_mit_gmx_file([ gsea_dir fname extname]);
  sets=match_with_affy(sets,gacc);
  [sets_lst,sets_idx]= get_gene_lists(sets);
  save([ gsea_dir fname '.mat'],'sets','sets_lst','sets_idx');
end

