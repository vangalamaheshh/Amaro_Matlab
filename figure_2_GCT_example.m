%make figure 2 example

SNPS=load_struct('~/Downloads/TCG-Testes_DFCI_7-Tumor-SM-4PDDH.snps.seg');
CPs=load_struct('~/Downloads/TCG-Testes_DFCI_7-Tumor-SM-4PDDH.CN.probes.seg');
SEG=load_struct('~/Downloads/TCG-Testes_DFCI_7-Tumor-SM-4PDDH.tsv');
SNPS.xstart=xhg19(chromosome2num_legacy(SNPS.chrom),str2double(SNPS.locstart));
CPs.xstart=xhg19(chromosome2num_legacy(CPs.contig),str2double(CPs.start));
SEG.xstart=xhg19(chromosome2num_legacy(SEG.Chromosome), str2double(SEG.Startbp));
SEG.xend=xhg19(chromosome2num_legacy(SEG.Chromosome),str2double(SEG.Endbp));
CPs.TCGTestes_DFCI_7TumorSM4PDDH=str2double(CPs.TCGTestes_DFCI_7TumorSM4PDDH);
SNPS.segmean=str2double(SNPS.segmean);
SEG.f=str2double(SEG.f);
flip_var=1;
figure()
hold on
for i=1:slength(SEG)
    ccp=reorder_struct(CPs,CPs.xstart>SEG.xstart(i)&CPs.xstart<SEG.xend(i));
    if flip_var==1
        color=[230/255,97/255,1/255];
        flip_var=0;
        plot(ccp.xstart(1:10:end),ccp.TCGTestes_DFCI_7TumorSM4PDDH(1:10:end),'b.')
    else
        flip_var=1;
        color=[94/255,60/255,153/255];
        plot(ccp.xstart(1:10:end),ccp.TCGTestes_DFCI_7TumorSM4PDDH(1:10:end),'r.')
    end
    
end
aL=num2chrom(1:23);
xL=xhg19(aL,zeros(size(aL)));
xM=round(conv(xL,[1 1]/2, 'same')); xM(end)=[];
set(gca,'xtick',xM,'xticklabel',aL,'xlim',[1 max(xL)])
for i=1:23
    line([xL(i),xL(i)],[0,10],'color',0.25*[1 1 1],'linestyle','--','LineWidth',.5)
end

figure()
hold on
for i=1:slength(SEG)
        csnp=reorder_struct(SNPS,SEG.xstart(i)<SNPS.xstart&SNPS.xstart<SEG.xend(i));
        if SEG.f(i) < .4 && ~isnan(SEG.f(i))
            plot(csnp.xstart(csnp.segmean<.5),csnp.segmean(csnp.segmean<.5),'.','Color',[230/255,97/255,1/255])
            plot(csnp.xstart(csnp.segmean>.5),csnp.segmean(csnp.segmean>.5),'.','Color',[94/255,60/255,153/255])

        else
            plot(csnp.xstart,csnp.segmean,'.','Color',[.5,.5,.5])
        end
end


aL=num2chrom(1:23);
xL=xhg19(aL,zeros(size(aL)));
xM=round(conv(xL,[1 1]/2, 'same')); xM(end)=[];
set(gca,'xtick',xM,'xticklabel',aL,'xlim',[1 max(xL)])
for i=1:23
    line([xL(i),xL(i)],[0,1],'color',0.25*[1 1 1],'linestyle','--','LineWidth',.5)
end
