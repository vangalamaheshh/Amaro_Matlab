%make whole genome bed file for aaron C. 

% Cytoband reference from here: http://genome.ucsc.edu/cgi-bin/hgTables
cytoband_ref=load_table('/Users/amaro/Downloads/FullGenomeCytoband.txt');
cytoband_ref=rmfield(cytoband_ref,'headline');
cytoband_ref=rmfield(cytoband_ref,'header');


chr.names=unique(cytoband_ref.chrom);
for i=1:length(chr.names)
    s=reorder_struct(cytoband_ref,ismember(cytoband_ref.chrom,chr.names{i}));
    chr.start(i,1)=1;
    chr.end(i,1)=max(s.chromEnd);
end
chr.names=chrom2num(chr.names);

chr=sort_struct(chr,'names');
chr=reorder_struct(chr,chr.names<23);

count=0;

    i=1;
    i
    bed.start=(chr.start(i):200:chr.end(i));
    bed.end=bed.start+100;
    bed.end=bed.end';
    bed.start=bed.start';
    bed.chr=zeros(length(bed.start),1,'single')+chr.names(i);
    BED=bed;

for i=2:slength(chr)
    i
    bed.start=(chr.start(i):200:chr.end(i));
    bed.end=bed.start+100;
    bed.end=bed.end';
    bed.start=bed.start';
    bed.chr=zeros(length(bed.start),1,'single')+chr.names(i);
    BED=mergeStruct(BED,bed);
    
end
BED=rmfield(BED,'N');
for i=1:slength(BED)
    BED.target{i,1}=strcat('T',num2str(i));
end
BED=reorderstructure(BED,'chr','start','end','target');
save_struct_noheader(BED,'/Users/amaro/Downloads/WGS_RECAPSEG.bed');