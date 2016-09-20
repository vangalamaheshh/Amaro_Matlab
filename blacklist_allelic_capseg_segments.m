function blacklist_allelic_capseg_segments(segfilein,segfileout)

seg=load_struct(segfilein);
seg.Chromosome=str2double(seg.Chromosome);
seg.Startbp=str2double(seg.Startbp);
seg.Endbp=str2double(seg.Endbp);


%chr9:43,484,541-68,928,174 centromere
blacklist.chr=[9];
blacklist.start=[43484541];
blacklist.end=[68928174];
%if any segment has 50% reciprocal overlap with a blacklisted segment
%remove it and merge the segments around it
b=zeros(slength(seg),1);


for i=1:slength(blacklist)
    k=find(blacklist.chr(i)==seg.Chromosome);
  for j=1:length(k)
              sm=k(j);
          % make sure that the start is less than the end 
    if (seg.Startbp(sm)<blacklist.end(i) && seg.Endbp(sm)>blacklist.end(i))
      s=max([seg.Startbp(sm);blacklist.start(i)]);
      e=min([seg.Endbp(sm);blacklist.end(i)]);
      
      if ((e-s)/(seg.Endbp(sm)-seg.Startbp(sm))) > .5 && ...
      ((e-s)/(blacklist.end(i)-blacklist.start(i))) > .5
         b(sm)=1;
               
      end
    end
   end

    
end
sego=reorder_struct(seg,~b);
save_struct(sego,segfileout);
end

function example
segfilein='/Users/amaro/Downloads/CLL-MDAC-0011-Tumor-SM-3VIEB.tsv';
blacklist_allelic_capseg_segments(segfilein,'/Users/amaro/Downloads/CLL-MDAC-0011-Tumor-SM-3VIEB.blacklist.tsv');
    
    
end