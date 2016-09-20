function gs=make_gene_set(hashes,lst,aa4hash)

for i=1:length(lst)
  gs(i).text=lst{i};
  
  % is in description
  disp(gs(i).text)
  [ord,s,res,h]=search_hashes(gs(i).text,hashes,aa4hash);
  if ~isempty(ord)
    [u,ui,uj]=unique_keepord(strvcat(aa4hash(ord,1)),'rows');
    [suj,suji]=sort(uj);
    gs(i).ord=ord(suji);
    gs(i).s=s(suji);
  else
    gs(i).ord=[];
    gs(i).s=[];
  end
end
