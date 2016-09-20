function C=expand_zone(C)

if mean(C.raw(:,1))>1
  disp('taking log of C.raw');
  C.raw(C.raw<0.1)=0.1;
  C.raw=log2(C.raw)-1;
end

C.zonechr=derunlength(C.zonechr_rl);
for i=1:size(C.dat,2)
  tmp=nan(size(C.zonechr,1),1);
  for j=1:max(C.zonechr_rl{i}(:,3));
    inzone=find(C.zonechr(:,i)==j);
    tmp(inzone)=median(C.raw(inzone,i),1);
  end
  C.zonechr(:,i)=tmp;
end

C.zonegen=derunlength(C.zonegen_rl);
for i=1:size(C.dat,2)
  tmp=nan(size(C.zonegen,1),1);
  for j=1:max(C.zonegen_rl{i}(:,3));
    inzone=find(C.zonegen(:,i)==j);
    tmp(inzone)=median(C.raw(inzone,i),1);
  end
  C.zonegen(:,i)=tmp;
end

