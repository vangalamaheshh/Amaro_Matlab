function gl=chrpos2genomic_location(chrn,r)

if length(chrn)>1
  gl=strcat(cellstr(repmat('chr',size(r,1),1)),...
            num2chromosome(chrn)',cellstr(repmat(':',size(r,1),1)),...
            cellstr(num2str(r(:,1),'%d')));
else
  gl=strcat(cellstr(repmat('chr',size(r,1),1)),...
            num2chromosome(chrn),cellstr(repmat(':',size(r,1),1)),...
            cellstr(num2str(r(:,1),'%d')));
end

if size(r,2)==2
  gl=strcat(gl,...
            cellstr(repmat('-',size(r,1),1)),cellstr(num2str(r(:,2),'%d')));
end

