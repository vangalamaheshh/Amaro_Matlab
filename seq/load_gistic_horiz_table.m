function T = load_gistic_horiz_table(fname)
% function for loading GISTIC peaks from
%      amp_genes.conf_95.*.txt
%   or del_genes.conf_95.*.txt

x = load_struct(fname);
if ~strcmp(x.cytoband{1},'q value') ||...
      ~strcmp(x.cytoband{2},'residual q value') ||...
      ~strcmp(x.cytoband{3},'wide peak boundaries') ||...
      ~strcmp(x.cytoband{4},'genes in wide peak')
  error('unexpected format');
end
f = fieldnames(x);

% also get first line as regular text to restore .'s in cytobands
tmp = fopen(fname);
h = fgetl(tmp);
fclose(tmp);
cb = split(h,char(9));

cb(1)=[]; f(1)=[]; % remove "cytoband" column
if strcmp('',cb(end))
  cb(end)=[]; f(end)=[];  % remove empty final column
end

T=[];
for i=1:length(cb)
  T.cytoband{i,1} = cb{i};
  T.q_value{i,1} = x.(f{i}){1};
  T.residual_q_value{i,1} = x.(f{i}){2};
  T.wide_peak_boundaries{i,1} = x.(f{i}){3};
  tmp = x.(f{i})(4:end);
  tmp(strcmp('',tmp))=[];
  T.genes_in_wide_peak{i,1} = decell(stringsplice(tmp',','));
end

