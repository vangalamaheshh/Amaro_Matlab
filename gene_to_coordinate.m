function [chr pos]=gene_to_coordinate(Hugo_Symbol,Gene_Table)

chr=chromosome2num_legacy(unique(Gene_Table.chrom(ismember(Gene_Table.name2,Hugo_Symbol))));
pos=str2double(unique(Gene_Table.txStart(ismember(Gene_Table.name2,Hugo_Symbol))));
pos=pos(1);


end



function test
Gene_Table=load_struct('/Users/amaro/Downloads/RefSeqGenes.txt');
Hugo_Symbol='TP53';

gene_to_coordinate(Hugo_Symbol,Gene_Table);


end