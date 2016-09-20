function [ output_args ] = absolute_allelic_cn_plot( seg1, GD, varargin )

if isempty(varargin)
    index=0;
else
    index = varargin{1};
end
if ~isfield(seg1,'xstart') && ~isfield(seg1,'xend')
    seg1.xstart=xhg19(seg1.Chromosome,seg1.Startbp); seg1.xend=xhg19(seg1.Chromosome,seg1.Endbp);
end

for i=1:slength(seg1) %merge and remove NaN segs
    if isnan(seg1.modala1(i))||isnan(seg1.modala2(i))
        if i<slength(seg1) & i > 1
            seg1.xend(i-1)=seg1.xend(i);
        end
    end
end

seg1=reorder_struct(seg1,~(isnan(seg1.modala1)|isnan(seg1.modala2)));

dx=.5;
%figure()
hold on
load('rb_colomap.mat')
ff={'modala1','modala2'};
if GD>1
    allele_thresh=2;
else
    allele_thresh=1;
end
for j=1:slength(seg1)
    for a=1:2 %2 alleles
        x1=[ seg1.xstart(j) seg1.xstart(j) seg1.xend(j) seg1.xend(j) seg1.xstart(j)];
        y1=[index+a-dx index+a+dx index+a+dx index+a-dx index+a-dx];
        
        if seg1.(ff{a})(j)>allele_thresh
            if seg1.(ff{a})(j)>(allele_thresh+3)
                c=[256,0,0];
            else
            c=rb_colormap(round((seg1.(ff{a})(j)/(allele_thresh+3))*128)+128,:);
            end
        elseif seg1.(ff{a})(j)==allele_thresh
            c=[256,256,256];
        elseif seg1.(ff{a})(j)<allele_thresh
            c=[0 0 256];
        end
        patch(x1,y1,c./256,'edgecolor','none')
    end
end
aL=num2chrom(1:23);
xL=xhg19(aL,zeros(size(aL)));
xM=round(conv(xL,[1 1]/2, 'same')); xM(end)=[];
set(gca,'xtick',xM,'xticklabel',aL,'xlim',[1 max(xL)])
for i=1:23
    line([xL(i),xL(i)],[.5,2.5],'color',0.25*[1 1 1],'linestyle','--','LineWidth',.5)
end
set(gca,'ytick',[])
%daspect([1440516642.50000,1,1]);
%pbaspect([1,0.182795698924731,0.182795698924731])
%print('~/Documents/MATLAB/test_18.eps','-depsc')
end

