%Make CN table for Matt

 FFPE_Seg=load_struct('/Users/amaro/Documents/BarrettsREvisionFigures/AggregatedFFPE_abs.seg');
 
 GenesForMatt={'EGFR';'ERBB2';'ERBB3';'FGFR2';'VEGFA';'KRAS';'ALK';'MET';'CCND1';'CCNE1';'CDK6';'MYB';'GATA4';'GATA6';'MYC';'PIK3CA';'CTNNB1';'IRS2'};
 Samples=unique(FFPE_Seg.sample);
 GeneTable=load_struct('/Users/amaro/Downloads/RefSeqGenes.txt');
 FFPE_Seg.xs=xhg19(chromosome2num_legacy(FFPE_Seg.Chromosome),str2double(FFPE_Seg.Startbp));
 FFPE_Seg.xe=xhg19(chromosome2num_legacy(FFPE_Seg.Chromosome),str2double(FFPE_Seg.Endbp));

 for i=1:length(GenesForMatt)
   l=find(ismember(GeneTable.name2,GenesForMatt{i}),1);
   Gstart=str2double(GeneTable.txStart(l));
   Gend=str2double(GeneTable.txEnd(l));
   Gchr=chromosome2num_legacy(GeneTable.chrom(l));
   Gstart=xhg19(Gchr,Gstart);
   Gend=xhg19(Gchr,Gend);
   OutPutTable.Gene{i,1}=GenesForMatt{i};
   for s=1:length(Samples)
       s_seg=reorder_struct(FFPE_Seg,ismember(FFPE_Seg.sample,Samples{s}));
       if isempty(find(s_seg.xe>Gstart&s_seg.xe>Gend&s_seg.xs>Gstart&s_seg.xs<Gend, 1))
           OutPutTable.(tr(Samples{s},'-','_'))(i,1)=max([str2double(s_seg.corrected_total_cn{find(s_seg.xe>Gstart,1,'first')});str2double(s_seg.corrected_total_cn{find(s_seg.xe<Gend,1,'first')})]);
       else
           OutPutTable.(tr(Samples{s},'-','_'))(i,1)=str2double(s_seg.corrected_total_cn{s_seg.xe>Gstart&s_seg.xe>Gend&s_seg.xs>Gstart&s_seg.xs<Gend});
       end
   end
   
    
    
end