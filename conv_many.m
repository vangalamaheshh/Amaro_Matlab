function c=conv_many(vecs)

if iscell(vecs)
  c=vecs{1};
  for i=2:length(vecs)
    c=conv(c,vecs{i});
  end
else
  [tmp,dc]=max(flipud(diff([vecs; zeros(1,size(vecs,2))],1,1)~=0),[],1);
  dc=size(vecs,1)+1-dc;
  
  c=vecs(1:dc(1),1);
  for i=2:size(vecs,2)
    c=conv(c,vecs(1:dc(i),i));
  end  
end
