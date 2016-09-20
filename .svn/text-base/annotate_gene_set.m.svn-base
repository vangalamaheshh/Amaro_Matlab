function annotate_gene_set(fname_in,fname_out,aa,hashes);

if strcmp(fname_in((end-2):end),'txt')
  pw_xls=read_dlm_file(fname_in,char(9));
  pw_xls=cat(1,pw_xls{:});
else
  [n,pw_xls,raw]=xlsread(fname_in);
  pw_xls=make_text_xlsread(raw);
end

if isempty(pw_xls{3,1}) % use second column
  tt=pw_xls(2:end,:);
  xx=cell(size(tt,1),1);
  for i=1:size(tt,1)
    xx{i}=sprintf(repmat('%s ; ',1,size(tt,2)),tt{i,:});
  end
  pw=make_gene_set(hashes,xx,[aa.symb aa.desc]);
else % use first
  pw=make_gene_set(hashes,pw_xls(2:end,1),[aa.symb aa.desc]);
end

if ~isfield(aa,'chip')
    warning('add chip field');
    aa.chip=cellstr(repmat('Chip',length(aa.symb),1));
end 

% FIX ME to work with any chip
f=fopen(fname_out,'w');
for i=1:length(pw)
  if isempty(pw(i).ord)
        fprintf(f,'%s\n',pw(i).text);
  else
    for j=1:min(length(pw(i).ord),10)
      if isfield(aa,'refgene')
        refgene=aa.refgene{pw(i).ord(j)};
        if isempty(refgene)
          rg='';
        else
          rg=refgene{1}{1};
        end
      else
        rg='';
      end
      % ds=aa.data{pw(i).ord(j),14};
      fprintf(f,'%s\t%d\t%s\t%s\t%s\t%s\t%s',...
              pw(i).text,j,'',aa.symb{pw(i).ord(j)},aa.desc{pw(i).ord(j)},...
              aa.chip{pw(i).ord(j)},aa.probeset{pw(i).ord(j)});
      if isfield(aa,'refgene')
        refgene=aa.refgene{pw(i).ord(j)};
        if isempty(refgene)
          rg='';
        else
          rg=refgene{1}{1};
        end
        fprintf(f,'\t%s',rg);
      end 
      fprintf(f,'\n');
    end
  end
end
fclose(f);
  
