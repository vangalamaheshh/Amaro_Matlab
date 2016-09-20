

ABS_seg_file=load_struct('/Users/amaro/Downloads/BarrettsCCFDp/Aggregate_Seg_file.txt');
sig_genes=load_struct('/Users/amaro/Downloads/genes.txt');
samples=load_struct('SampleTableBarretts_TSPs.txt');



Gene_Table=load_struct('/Users/amaro/Downloads/RefSeqGenes.txt');
absolute_table=load_struct('/Users/amaro/Downloads/BarrettsCCFDp/5.21.2014_redo.ATW.ABSOLUTE.table.txt');

ABS_maf=load_struct('/Users/amaro/Downloads/BarrettsCCFDp/Bmaffilt.txt');

ABS_seg_file.Chromosome=str2double(ABS_seg_file.Chromosome);
ABS_seg_file.Endbp=str2double(ABS_seg_file.Endbp);
ABS_seg_file.Startbp=str2double(ABS_seg_file.Startbp);
ABS_seg_file.total_HZ=str2double(ABS_seg_file.total_HZ);
ABS_seg_file.HZ=str2double(ABS_seg_file.HZ);
ABS_seg_file.SC_HZ=str2double(ABS_seg_file.SC_HZ);
ABS_seg_file.corrected_total_cn=str2double(ABS_seg_file.corrected_total_cn);
%
absolute_table.ploidy=str2double(absolute_table.ploidy);

for i=1:slength(samples)
k=find(ismember(absolute_table.array,samples.pair_id{i}));
samples.GD{i,1}=absolute_table.Genomedoublings{k};
end
for v=1:slength(sig_genes)
for s=1:slength(samples)
    [abs_cn,HZ]=get_copy_number_of_gene(samples.pair_id{s},ABS_seg_file,sig_genes.gene{v},Gene_Table);
    ploidy=absolute_table.ploidy(ismember(absolute_table.array,samples.pair_id{s}));
    k=ismember(ABS_maf.sample,samples.names{s});
if ismember(sig_genes.gene{v},ABS_maf.Hugo_Symbol(k))
    samples.(sig_genes.gene{v}){s,1}='Mut';
elseif HZ==1 && v <8
    samples.(sig_genes.gene{v}){s,1}='HZ';
elseif (abs_cn/ploidy) >2 && v >= 8
    samples.(sig_genes.gene{v}){s,1}='AMP';

else
    samples.(sig_genes.gene{v}){s,1}='None';
end
end

sig_genes.gene{v}
end
save_struct(samples,'Barretts_Samples_Events.txt');