function cmat=chi2_box_v_sbox(box,sbox)

for i=1:size(box,1)
  for j=1:size(sbox,1)
    [t,c2,p]=crosstab(box(i,:),sbox(j,:));
    cmat(i,j)=p;
  end
  i
end

