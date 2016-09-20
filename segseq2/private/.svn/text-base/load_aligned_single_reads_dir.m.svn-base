function READS=load_aligned_single_reads_dir( dirname, qualCutoff, sampleName )

matfile = [ sampleName '_aligned_reads_qual' num2str(qualCutoff) '.mat' ];
chrList = 1:23;

READS.chr = [];
READS.pos1 = [];
READS.lane = [];

for c=1:length(chrList)
    chr=chrList(c);
    filename = [ dirname '/chr' num2str(chr) '.txt' ]
    R = load_aligned_single_reads( filename, qualCutoff );

    READS.chr = [ READS.chr; R.chr ];
    READS.pos1 = [ READS.pos1; R.pos1 ];
    READS.lane = [ READS.lane; R.lane ];
end

save(matfile,'READS','-v7.3');
