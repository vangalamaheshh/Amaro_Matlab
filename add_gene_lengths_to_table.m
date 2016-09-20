%add lengths to gene table

Gene_Table=load_struct('/Users/amaro/Downloads/RefSeqGenes.txt');
for i=1:slength(Gene_Table)

strs=split(Gene_Table.exonStarts{i},',');
strs{1}=strs{1}(2:end);

ens=split(Gene_Table.exonEnds{i},',');
ens{1}=ens{1}(2:end);


Gene_Table.Length(i,1)=sum(str2double({ens{1:end-1}})-str2double({strs{1:end-1}}));

end