function annotate_gene_set(fname_in,fname_out,aa,hashes);

[n,pw_xls]=xlsread(fname_in);

if isempty(pw_xls{2,1}) % use second column
  pw=make_gene_set(hashes,pw_xls(2:end,2),[aa.symb aa.desc]);
else % use first
  pw=make_gene_set(hashes,pw_xls(2:end,1),[aa.symb aa.desc]);
end

% FIX ME to work with any chip
f=fopen(fname_out,'w');
for i=1:length(pw)
  for j=1:min(length(pw(i).ord),10)
    refgene=aa.refgene{pw(i).ord(j)};
    if isempty(refgene)
      rg='';
    else
      rg=refgene{1}{1};
    end
    ds=aa.data{pw(i).ord(j),14};
    fprintf(f,'%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n',...
            pw(i).text,j,'',aa.probeset{pw(i).ord(j)},...
            aa.symb{pw(i).ord(j)},rg,ds,aa.desc{pw(i).ord(j)});
  end
end
fclose(f);
