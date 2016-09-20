function L = load_concordance_file(fname)

demand_file(fname)
L = load_lines(fname);
idx = find(~cellfun('isempty',L)&~strncmp(L,'#',1),1);
if isempty(idx)
  L=[];
else
  L = load_struct_specify_numeric_cols(fname,[],idx,true);
  L = make_numeric(L,{'reference','non_reference','pct_concordance'});
end
