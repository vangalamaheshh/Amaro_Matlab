function cm=make_colormap(origcm,n)

for i=1:3
  cm(:,i)=resample(origcm(:,i),n,64);
end

cm(cm<0)=0;
cm(cm>1)=1;

