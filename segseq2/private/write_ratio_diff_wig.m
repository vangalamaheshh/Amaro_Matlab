function write_ratio_diff_wig( fname, R, CHR, POS, stepIdx )

fid = fopen( fname, 'w' );

fprintf(fid,'track type=wiggle_0\n');
chrList = unique(CHR);
for c=1:length(chrList)
    idxChr = find(CHR==chrList(c));
    chrStart = POS(idxChr);
    chrDiff = R(idxChr);

    [sortPos,idxOrder] = sort(chrStart);
    sortDiff = chrDiff(idxOrder);

    if chrList(c) == 23
        fprintf(fid,'variableStep chrom=chrX\n');
    else
        fprintf(fid,'variableStep chrom=chr%d\n', chrList(c));
    end
    
    idxToPrint=1:stepIdx:length(sortPos);
    for i=1:length(idxToPrint)
	fmt=['%d\t%.3f\n'];
        fprintf(fid,fmt,sortPos(idxToPrint(i)),sortDiff(idxToPrint(i)));
    end
end
fclose(fid);
