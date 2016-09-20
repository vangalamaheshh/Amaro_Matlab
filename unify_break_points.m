%Tau/2-delta/2=mu minor
%delta=(mu_major-mu_minor)/2

%weighting segments w = length/sigmadelta^2
%2sigma^2=sigmaDelta^2
%every segment has delta

%seperate function take the union of the breakpoints and segment the data equally
% then the fit
%ex fit type
% g = fittype('a*x')
% linfit=fit(call.,call.normal_f,g,'Weights',w)
%
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

function test
A_segs=load_struct('/xchip/cga_home/amaro/CLL/Sample_QC/LOH/a_segments.tsv');
for i=1:slength(A_segs)
TumorSeg=load_table(A_segs.alleliccapseg_tsv{i});
NormalSeg=load_table(A_segs.allelic_capseg_seg_normal{i});


TumorSeg.delta=TumorSeg.mu_major-TumorSeg.mu_minor;
NormalSeg.delta=NormalSeg.mu_major-NormalSeg.mu_minor;

TumorSeg.sigma=2*(TumorSeg.sigma_major.^2);
NormalSeg.sigma=2*(NormalSeg.sigma_major.^2);
[Tseg Nseg]=unify_break_points(TumorSeg,NormalSeg);

%weighting segments w = length/sigmadelta^2
Tseg.length=Tseg.End-Tseg.Start;
Nseg.length=Nseg.End-Nseg.Start;
Tseg=reorder_struct(Tseg,~(Tseg.Chromosome==23));
Nseg=reorder_struct(Nseg,~(Nseg.Chromosome==23));
Tseg=reorder_struct(Tseg,(Tseg.length>1*10^6));
Nseg=reorder_struct(Nseg,(Nseg.length>1*10^6));

W=Tseg.length./(Tseg.sigma.^2+Nseg.sigma.^2);
Tseg=reorder_struct(Tseg,~isnan(W));
Nseg=reorder_struct(Nseg,~isnan(W));

W=W(~isnan(W));
%ex fit type
 g = fittype('a*x');
 linfit=fit(Tseg.delta,Nseg.delta,g,'Weights',W);
 hold on
 BubblePlot(Tseg.delta,Nseg.delta,Tseg.length,[0 0 1],300)
 hline=refline(linfit.a,0);
 set(hline,'Color',[0 0 0]);
  title(strcat(A_segs.pair_id{i},'_',num2str(linfit.a)));
saveas(gcf,strcat(A_segs.pair_id{i},'.fig'));
A_segs.LOH_estimate(i,1)=linfit.a;
close(gcf)


end




end
