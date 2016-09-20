function [] = writeFileWithSingleNumber(filename, value)
fid = fopen(filename, 'w');
fprintf(fid, '%d', value);
fclose(fid);