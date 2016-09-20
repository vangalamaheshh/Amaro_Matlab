function [map] = make_sample_table_ESO(ABS_maf,ABS_seg_file,sig_genes,absolute_table,map,old_sample_table,Gene_Table)
map=reorder_struct(map,ismember(map.pair_id,old_sample_table.pair_id)); %perserve blacklisting
%headers needed:
%names	pair_id	gistic_id	order	tissue	relatedness	Individual	GD	EGFR	ERBB2	ERBB3	FGFR2	VEGFA	KRAS	ALK	MET	CCND1	CCNE1	CDK6	MYB	GATA4	GATA6	MYC	CTNNB1
%genes should contain AMP/Mut/None
map.GD=zeros(slength(map),1);


ABS_seg_file.Chromosome=chromosome2num_legacy(ABS_seg_file.Chromosome);
ABS_seg_file.Endbp=str2double(ABS_seg_file.Endbp);
ABS_seg_file.Startbp=str2double(ABS_seg_file.Startbp);
ABS_seg_file.total_HZ=str2double(ABS_seg_file.total_HZ);
ABS_seg_file.HZ=str2double(ABS_seg_file.HZ);
ABS_seg_file.SC_HZ=str2double(ABS_seg_file.SC_HZ);
ABS_seg_file.corrected_total_cn=str2double(ABS_seg_file.corrected_total_cn);


[i m]=ismember(old_sample_table.pair_id,map.pair_id);
map.tissue(m(m>0),1)=old_sample_table.tissue(i);
map.relatedness(m(m>0),1)=old_sample_table.relatedness(i);
map.GD=zeros(slength(map),1);

for i=1:slength(sig_genes)
    map.(sig_genes.gene{i})=repmat({'None'},slength(map),1);
end

ABS_maf=reorder_struct(ABS_maf,ismember(ABS_maf.i_judgement,'KEEP'));
ABS_maf=reorder_struct(ABS_maf,ismember(ABS_maf.Variant_Classification,'Missense_Mutation')|ismember(ABS_maf.Variant_Classification,'Nonsense_Mutation')|...
ismember(ABS_maf.Variant_Classification,'Splice_Site')|ismember(ABS_maf.Variant_Classification,'Frame_Shift_Ins')|ismember(ABS_maf.Variant_Classification,'Frame_Shift_Del')|ismember(ABS_maf.Variant_Classification,'In_Frame_Del')|ismember(ABS_maf.Variant_Classification,'In_Frame_Ins'));
absolute_table.Genomedoublings=str2double(absolute_table.Genomedoublings);     

for i=1:slength(map)
    k=find(ismember(absolute_table.sample,map.pair_id{i}));
    map.GD(i)=absolute_table.Genomedoublings(k)>0;
end

for i=1:slength(sig_genes)
    for j=1:slength(map)
       [abs_cn,HZ]=get_copy_number_of_gene(map.pair_id{j},ABS_seg_file,sig_genes.gene{i},Gene_Table);
       map.pair_id{j}
       ploidy=str2double(absolute_table.ploidy(ismember(absolute_table.sample,map.pair_id{j})));
        k=ismember(ABS_maf.sample,map.case_sample{j});
         if ismember(sig_genes.gene{i},ABS_maf.Hugo_Symbol(k))
             map.(sig_genes.gene{i}){j}='Mut';
         end
         if HZ==1
             map.(sig_genes.gene{i}){j}='DEL';
         end
         
    end
end



end

function test
ABS_maf=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/TSPsMaf.txt');
sig_genes=load_struct('/Users/amaro/Downloads/BarrettsCCFDp/genes_for_ESO_TSPs.txt');
absolute_table=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/BarrettsABSOLUTEtable.txt');

map=load_struct('/Users/amaro/Downloads/BarrettsCCFDp/samplemapping.txt');
ABS_seg_file=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/Aggregate_Seg_file.txt');
old_sample_table=load_struct('Barretts_Samples_Events_onc_new_order.txt');
Gene_Table=load_struct('/Users/amaro/Downloads/RefSeqGenes.txt');
map=make_sample_table_ESO(ABS_maf,ABS_seg_file,sig_genes,absolute_table,map,old_sample_table,Gene_Table)
end