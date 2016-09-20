function rg=mark_futreal_genes(rg)

futreal=read_dlm_file('~gadgetz/projects/snp/data/Futreal/Table_1_full_2006-02-16.txt');
for i=2:length(futreal)
  fsymb{i-1}=futreal{i}{1};
  flocusid(i-1)=str2num(futreal{i}{3});
end

for i=1:length(rg)
  rg(i).symbol=rg(i).symb;
end

%  [Mt,m1,m2]=match_string_sets_hash(fsymb,{rg.symb});
%  missing1=setdiff(1:length(fsymb),unique(m1));
[Nt,n1,n2]=match_string_sets_hash(cellstr(num2str(flocusid')),cellstr(num2str([rg.locus_id]')));
missing1a=setdiff(1:length(fsymb),unique(n1));
un2=unique(n2);
for i=1:length(un2)
  if rg(un2(i)).symbol(end)~='*'
    rg(un2(i)).symbol=[ rg(un2(i)).symbol '**'];
  end
end
