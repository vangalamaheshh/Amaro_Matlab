function [uSeg1 uSeg2]=unify_break_points(seg1,seg2)

seg=mergeStruct(seg1,seg2);

xs1=xhg19(seg1.Chromosome,seg1.Start_bp);
xe1=xhg19(seg1.Chromosome,seg1.End_bp);

xs2=xhg19(seg2.Chromosome,seg2.Start_bp);
xe2=xhg19(seg2.Chromosome,seg2.End_bp);



x1=xhg19(seg.Chromosome,seg.Start_bp);
x2=xhg19(seg.Chromosome,seg.End_bp);

X=unique(sort([x1;x2]));
Starts=X(1:(end-1));
Ends=X(2:end)-1;

uSeg1.Start=Starts;
uSeg1.End=Ends;
uSeg1.Med=median(Starts,Ends);


uSeg2.Start=Starts;
uSeg2.End=Ends;
uSeg2.Med=median(Starts,Ends);

for i=1:length(Starts)
    e=find(uSeg1.Med(i)<=xe1,1,'First');
    s=find(uSeg1.Med(i)>=xs1,1,'Last');
    if ~s==e
        uSeg1.tau(i)=NaN;
        
    else
        uSeg1.Chromosome(i,1)=seg1.Chromosome(s);
        uSeg1.tau(i,1)=seg1.tau(s);
        uSeg1.delta(i,1)=seg1.delta(s);
        uSeg1.sigma(i,1)=seg1.sigma(s);
    end
    
    e=find(uSeg2.Med(i)<=xe2,1,'First');
    s=find(uSeg2.Med(i)>=xs2,1,'Last');
    if ~s==e
        uSeg2.tau(i)=NaN;
        
    else
        uSeg2.Chromosome(i,1)=seg2.Chromosome(s);
        uSeg2.tau(i,1)=seg2.tau(s);
        uSeg2.delta(i,1)=seg2.delta(s);
        uSeg2.sigma(i,1)=seg2.sigma(s);
    end
    
    
    
    
end

uSeg1=reorder_struct(uSeg1,~isnan(uSeg1.tau));
    
uSeg2=reorder_struct(uSeg2,~isnan(uSeg2.tau));

end