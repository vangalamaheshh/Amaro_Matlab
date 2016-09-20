function make_gmt_from_symb(gmt_fname,symb)

f=fopen(gmt_fname,'w');
if ~iscell(symb)
  symb={symb};
end

for i=1:length(symb)
  s=read_dlm_file(symb{i},char(9));
  s=cat(1,s{:});
  symbol_col=find(~cellfun('isempty',regexp(lower(s(1,:)),'symbol')));
  if length(symbol_col)~=1
    error('No unique symbol column');
  end
  set_name=regexprep(symb{i},'\.symb\.txt','');
  ref_name=s{3,1};
  if isempty(ref_name)
    ref_name='BLACK';
  end
  fprintf(f,'%s\t%s',set_name,ref_name);
  for j=2:size(s,1)
    if ~strcmp(s{j,symbol_col},'---')
      fprintf(f,'\t%s',upper(s{j,symbol_col}));
    end
  end
  fprintf(f,'\n');
end
fclose(f);
