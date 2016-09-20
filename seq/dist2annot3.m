function dist=dist2annot3(bkpts,annotpos,chrinds)

n=length(bkpts.chr1);
dist=zeros(n,2);
for i=1:n
    dist(i,1)=min(abs(bkpts.pos1(i)-empty2nan(annotpos(chrinds{bkpts.chr1(i)}))));
    dist(i,2)=min(abs(bkpts.pos2(i)-empty2nan(annotpos(chrinds{bkpts.chr2(i)}))));    
end
