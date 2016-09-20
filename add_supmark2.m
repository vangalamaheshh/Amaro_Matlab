function D=add_supmark2(D,c)

for i=1:size(D.supdat,1)
  D.supmark(i).colormap=c(1:floor((size(c,1)-1)/(length(unique(D.supdat(i,:)))-1)):size(c,1),:);
  D.supmark(i).height=1;
  D.supmark(i).patchwidth=1;      
  D.supmark(i).linewidth=0;
end
