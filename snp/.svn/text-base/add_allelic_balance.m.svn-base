function C=add_allelic_balance(C,L)

for i=1:length(C.peaks)
  cab=zeros(1,length(C.peaks{i}));
  for j=1:length(C.peaks{i})
    in_pk=find(C.level(:,i)==j);
    if ~isempty(in_pk)
      cab(j)=sum(L.dat(in_pk,i)==0,1)./length(in_pk);
    end
  end
  C.peaks_allele_balance{i}=cab;
end

