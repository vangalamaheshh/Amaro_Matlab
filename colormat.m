function colormat(dat,convtab)

cmat=zeros(size(dat,1),size(dat,2),3);
for i=1:size(dat,1)
  for j=1:size(dat,2)
    if iscell(dat)
      cmat(i,j,:)=convert_enum(dat{i,j},convtab);
    else
      cmat(i,j,:)=convert_enum(dat(i,j),convtab);
    end      
  end
end

image(cmat);

