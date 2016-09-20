function C=fix_C(C)

C=shrink_D(C);

if mean(C.raw(:,1))>1
  disp('taking log of C.raw');
  C.raw(C.raw<0.1)=0.1;
  C.raw=log2(C.raw)-1;
end

disp('fixing rl segments');
for i=1:size(C.dat,2)
  cur=C.cbs_rl{i};
  for j=1:size(cur,1)
    if cur(j,2)<=size(C.raw,1)
      cur(j,3)=median(C.raw(cur(j,1):cur(j,2),i),1);
    else
      cur(j,3)=NaN;
    end
  end
  cur(isnan(cur(:,3)),:)=[];
  C.cbs_rl{i}=cur;
  fprintf(1,'.');
end
disp('Done');

disp('fixing dat');
C.dat=derunlength(C.cbs_rl);

disp('calculate medians');
C.medians=median(C.dat(C.chrn~=23,:),1);

disp('sub. median from dat');
C.dat=C.dat-repmat(C.medians,size(C.dat,1),1);
disp('sub. median from raw');
C.raw=C.raw-repmat(C.medians,size(C.raw,1),1);

disp('sub. median from rl');
for i=1:size(C.dat,2)
  C.cbs_rl{i}(:,3)=C.cbs_rl{i}(:,3)-C.medians(i);
end
