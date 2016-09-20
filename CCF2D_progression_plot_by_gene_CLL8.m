sigABS=load_table('~/Documents/SuppTable9WithValidationCountsAdditionalSamples.tsv');
sig_genes=load_struct('/Volumes/xchip_cga_home/amaro//CLL/StilgenbauerMafFreeze2.0/ICGC_StilgenbauerForComut1.16/ICGC_Stilgenbauer.aggregate.NOIGV_WITHTTKandNOTCH1_PoNCut4_WithSalvage.sig_genes_.txt');
M=load_struct('~/Documents/FullSetCLL8_CCFs.txt');
M=reorder_struct(M,~cellfun(@isempty,strfind(M.Tumor_Sample_Barcode,'GCLL')));
XM=reorder_struct(M,ismember(M.Chromosome,'X'));
pair_map=load_struct('/Volumes/xchip_cga_home/amaro/CLL/Stilgenbauer.DFCI.ICGC.630.pairs.tsv');
for i=1:slength(XM)
k=find(ismember(pair_map.case_sample,XM.Tumor_Sample_Barcode{i}));
XM.pair_id{i,1}=pair_map.pair_id{k};
end



external_ids=load_struct('/Volumes/xchip_cga_home/amaro/CLL/ABSOLUTE/full_list_of_external_samples.tsv');
external_ids=reorder_struct(external_ids,~ismember(external_ids.external_id_capture,''));
cheatsheet=load_struct('/Users/amaro/Downloads/c.txt');
for i=1:slength(sigABS)
sigABS.indi{i,1}=sigABS.sample1{i}(5:13);
end
sigABS=reorder_struct(sigABS,ismember(sigABS.indi,cheatsheet.sample));


for i=1:slength(external_ids)
external_ids.TP{i,1}=external_ids.external_id_capture{i}(end-1:end);
external_ids.sm_id{i,1}=external_ids.sample_id{i}(end-7:end);
end

%sigABS=reorder_struct(sigABS,ismember(sigABS.locus,sig_genes.gene)); %keep
%maf intact for X chrom 

for i=1:slength(sigABS)
for j=1:slength(external_ids)
if ~isempty(strfind(sigABS.sample2{i},external_ids.sm_id{j}))
sigABS.TP2{i,1}=external_ids.TP{j};
end
end
end
sigABS=reorder_struct(sigABS,ismember(sigABS.TP2,'02'));


for i=1:slength(pair_map)
    XM.pair_id(ismember(XM.Tumor_Sample_Barcode,pair_map.case_sample{i}),1)=pair_map.pair_id(i);
end

XM=reorder_struct(XM,ismember(XM.pair_id,sigABS.sample1)|ismember(XM.pair_id,sigABS.sample2));
for i=1:slength(XM)
    XM.individual{i,1}=XM.pair_id{i}(5:13);
    k=find(ismember(external_ids.sample_id,XM.Tumor_Sample_Barcode{i}));
    XM.TP{i,1}=external_ids.TP{k};

end
XM.ccf_median_hack=str2double(XM.ccf_median_hack);
XM.t_alt_count=str2double(XM.t_alt_count);
patients_with_x=unique(XM.individual);
counter=1;
for i=1:length(patients_with_x)
p_sigABS=reorder_struct(sigABS,ismember(sigABS.indi,patients_with_x{i}));
cluster_centers=unique([p_sigABS.CCF1,p_sigABS.CCF2],'rows','stable');
highcent=unique([p_sigABS.CCF1high,p_sigABS.CCF2high],'rows','stable');
lowcent=unique([p_sigABS.CCF1low,p_sigABS.CCF2low],'rows','stable');

p_XM=reorder_struct(XM,ismember(XM.individual,patients_with_x{i}));

sites_to_assign=unique(p_XM.Start_position);
for j=1:length(sites_to_assign)
    TP1_ccf=max([p_XM.ccf_median_hack((ismember(p_XM.TP,'01')&ismember(p_XM.Start_position,sites_to_assign{j}))) 0]);
    TP2_ccf=max([p_XM.ccf_median_hack((ismember(p_XM.TP,'02')&ismember(p_XM.Start_position,sites_to_assign{j}))) 0]);
    d_m=(pdist([TP1_ccf,TP2_ccf;cluster_centers],'euclidean'));
    [trsh, indx]=min(d_m(1:length(cluster_centers)));
    xABS.CCF1(counter,1)=cluster_centers(indx,1);
    xABS.CCF2(counter,1)=cluster_centers(indx,2);
    xABS.CCF1high(counter,1)=highcent(indx,1);
    xABS.CCF1low(counter,1)=lowcent(indx,1);
    xABS.CCF2high(counter,1)=highcent(indx,2);
    xABS.CCF2low(counter,1)=lowcent(indx,2);
    xABS.locus(counter,1)=unique(p_XM.Hugo_Symbol(ismember(p_XM.Start_position,sites_to_assign{j})));
    xABS.Chr{counter,1}='X';
    xABS.pos{counter,1}=sites_to_assign{j};
    xABS.Reference_Allele(counter,1)=unique(p_XM.Reference_Allele(ismember(p_XM.Start_position,sites_to_assign{j})));
    xABS.Mutant_Allele(counter,1)=unique(p_XM.Tumor_Seq_Allele2(ismember(p_XM.Start_position,sites_to_assign{j})));
    xABS.Protein_Change(counter,1)=unique(p_XM.Protein_Change(ismember(p_XM.Start_position,sites_to_assign{j})));
    xABS.Variant_Type(counter,1)=unique(p_XM.Variant_Type(ismember(p_XM.Start_position,sites_to_assign{j})));
    xABS.Variant_Classification(counter,1)=unique(p_XM.Variant_Classification(ismember(p_XM.Start_position,sites_to_assign{j})));
    if ~isempty(unique(p_XM.pair_id((ismember(p_XM.TP,'01')))))
    xABS.sample1(counter,1)=unique(p_XM.pair_id((ismember(p_XM.TP,'01'))));
    else
        xABS.sample1{counter,1}='NaN';
    end
    if ~isempty(unique(p_XM.pair_id((ismember(p_XM.TP,'02')))))
    xABS.sample2(counter,1)=unique(p_XM.pair_id((ismember(p_XM.TP,'02'))));
    else
        xABS.sample2{counter,1}='NaN';
    end

    counter=counter+1;
    
    
end

end

xABSfields=fieldnames(xABS);
sigABSfields=fieldnames(sigABS);
sigABSfields(ismember(sigABSfields,xABSfields));

sigABS=rmfield(sigABS,sigABSfields(~ismember(sigABSfields,xABSfields)));
sigABS=mergeStruct(sigABS,xABS);
sigABS=rmfield(sigABS,'N');
sigABS=reorder_struct(sigABS,ismember(sigABS.locus,sig_genes.gene));
genes=unique(sigABS.locus);

for s=1:length(genes)
   g_maf=reorder_struct(sigABS,ismember(sigABS.locus,genes{s})); 
   
    
 g_maf.CCF1low=g_maf.CCF1-g_maf.CCF1low;
 g_maf.CCF2low=g_maf.CCF2-g_maf.CCF2low;
 g_maf.CCF1high=g_maf.CCF1high-g_maf.CCF1;
 g_maf.CCF2high=g_maf.CCF2high-g_maf.CCF2;
tp1x=.8:.4/slength(g_maf):1.2;
tp2x=1.8:.4/slength(g_maf):2.2;
g_maf=sortstruct(g_maf,'CCF1',-1);
tp1x=[.8:.4/(slength(g_maf)-1):1.2]';
g_maf=sortstruct(g_maf,'CCF2',-1);
tp2x=[1.8:.4/(slength(g_maf)-1):2.2]';

num_up=sum((g_maf.CCF1+g_maf.CCF1high)<(g_maf.CCF2-g_maf.CCF2low));
num_down=sum((g_maf.CCF1-g_maf.CCF1low)>(g_maf.CCF2+g_maf.CCF2high));
total=slength(g_maf);
if slength(g_maf)<3
    p_up=1;
    p_down=1;
else
p_up=1-binocdf(num_up,total,.5);
p_down=1-binocdf(num_down,total,.5);


p.shifted(s,1)=1-binocdf(num_up+num_down,total,.5);
p.unshifted(s,1)=1-binocdf(total-(num_up+num_down),total,.5);
p.gene{s,1}=genes{s};
p.num_up(s,1)=num_up;
p.num_down(s,1)=num_down;
p.stable(s,1)=total-(num_up+num_down);


end

if num_up+num_down<4
    p.shifted(s,1)=1-binocdf(num_up+num_down,total,.5);
p.unshifted(s,1)=1-binocdf(total-(num_up+num_down),total,.5);
    p_up_l=1;
    p_down_l=1;
    p.gene{s,1}=genes{s};
p.num_up(s,1)=num_up;
p.num_down(s,1)=num_down;
p.stable(s,1)=total-(num_up+num_down);

else
p_up_l=1-binocdf(num_up,num_up+num_down,.5);
p_down_l=1-binocdf(num_down,num_up+num_down,.5);
end



   hold on   
for i=1:slength(g_maf)
       if (g_maf.CCF1(i)+g_maf.CCF1high(i))<(g_maf.CCF2(i)-g_maf.CCF2low(i))
        errorbar([tp1x(i);tp2x(i)],[g_maf.CCF1(i);...
            g_maf.CCF2(i)],[g_maf.CCF1low(i);g_maf.CCF2low(i)],...
        [g_maf.CCF1high(i);g_maf.CCF2high(i)],'r')
       elseif (g_maf.CCF1(i)-g_maf.CCF1low(i))>(g_maf.CCF2(i)+g_maf.CCF2high(i))
                             errorbar([tp1x(i);tp2x(i)],[g_maf.CCF1(i);...
            g_maf.CCF2(i)],[g_maf.CCF1low(i);g_maf.CCF2low(i)],...
        [g_maf.CCF1high(i);g_maf.CCF2high(i)],'b')
       else
                   errorbar([tp1x(i);tp2x(i)],[g_maf.CCF1(i);...
            g_maf.CCF2(i)],[g_maf.CCF1low(i);g_maf.CCF2low(i)],...
        [g_maf.CCF1high(i);g_maf.CCF2high(i)],'Color',[.5 .5 .5])

       end
               plot([tp1x(i)],g_maf.CCF1(i),'k.')
               plot([tp2x(i)],g_maf.CCF2(i),'k.')
               text(2.3, .7, strcat('p_total=',num2str(min([p_up p_down]))));
               text(2.3, .5, strcat('p_up_down=',num2str(min([p_up_l p_down_l]))));

end
   
xlim([0.5 2.5])
ylim([0 1]);
title(strcat(genes{s}),'FontSize',28)
ylabel('Cancer Cell Fraction','FontSize',20)
set(gca,'YTick',[0;.5;1],'XTick',[1;2],'XTickLabel',{'TP1';'TP2'})
print(gcf,'-depsc',strcat('/Users/amaro/Documents/RevisionFiguresCLL8/MutationProgressionPlots/',genes{s},'.eps'));
close all


end

