function [ch,ah]=plot_chr_arm_boundaries(C,cyto)
  
if ~isfield(C,'armn')
  C=add_cyto(C,cyto);
end

chrpos=find(diff(C.chrn));
armpos=find(diff(C.armn)>0);

ax=axis;
for i=1:length(chrpos)
  ch(i)=line([ chrpos(i) chrpos(i)],ax(3:4),'Color','r');
end
for i=1:length(armpos)
  ah(i)=line([ armpos(i) armpos(i)],ax(3:4),'Color','g');
end
