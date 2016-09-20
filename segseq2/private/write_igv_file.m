function write_seg_file( fname, SEG, sampleName )

fid = fopen( fname, 'w' );

fmt=['Sample\tChromosome\tStart\tEnd\tLog2 copy ratio\n'];
fprintf(fid,fmt);
%fmt=['%s\t%d\t%.0f\t%.0f\t%.2e\t%.4f\n'];
fmt=['%s\t%d\t%.0f\t%.0f\t%.4f\n'];

for i=1:length(SEG.chr)
	fprintf(fid, fmt, sampleName, SEG.chr(i), SEG.left(i), SEG.right(i), log2(SEG.ratios(i)) );
end
fclose(fid);
