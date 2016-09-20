function [cp_num HZ]=get_copy_number_of_gene(sample_id,ABS_seg_file,gene,Gene_Table)

[chr pos]=gene_to_coordinate(gene,Gene_Table);
k=ismember(ABS_seg_file.sample,sample_id)&ABS_seg_file.Chromosome==chr;
l=find(k);
for i=1:length(l)
    
        if ABS_seg_file.Startbp(l(i))<=pos && ABS_seg_file.Endbp(l(i))>=pos
            cp_num=ABS_seg_file.corrected_total_cn(l(i));
            HZ=max([ABS_seg_file.total_HZ(l(i)),ABS_seg_file.HZ(l(i)),ABS_seg_file.SC_HZ(l(i))]);
            return 
        end
        if ABS_seg_file.Endbp(l(i))>=pos
            cp_num=ABS_seg_file.corrected_total_cn(l(i));
            HZ=max([ABS_seg_file.total_HZ(l(i)),ABS_seg_file.HZ(l(i)),ABS_seg_file.SC_HZ(l(i))]);

            return 
        
        end




end
end

function test
Gene_Table=load_struct('/Users/amaro/Downloads/RefSeqGenes.txt');
ABS_seg_file=load_struct('/Users/amaro/Downloads/BarrettsCCFDp/Aggregate_Seg_file.txt');

sample_id='UMBEER_101-TP-NWGA-SM-2TWM-SM-2TSI';
gene='ERBB2'
get_copy_number_of_gene('UMBEER_101-TP-NWGA-SM-2TWM-SM-2TSI',ABS_seg_file,'ERBB2',Gene_Table)

end