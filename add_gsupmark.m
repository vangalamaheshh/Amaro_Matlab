function D=add_gsupmark(D,c)

for i=1:size(D.gsupdat,1)
  [st,sn]=break_sup_names(deblank(D.gsupacc(i,:)));
  m=ceil(max(D.gsupdat(i,:),[],2));
  if ~isempty(sn) | m>2
    c1=repmat(c,max(length(sn),m),1);
    D.gsupmark(i).colormap=c1(1:max(length(sn),m),:);
  else
    D.gsupmark(i).colormap=[0 1 0; 1 0 0];
  end
  D.gsupmark(i).height=1;
end
