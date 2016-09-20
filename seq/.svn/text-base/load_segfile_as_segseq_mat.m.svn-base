function S = load_segfile_as_segseq_mat(fname)

z = load_struct(fname);
S.RATIOS.chr = convert_chr(z.chromosome);
S.RATIOS.windows = round((str2double(z.start) + str2double(z.end))/2);
S.RATIOS.ratios = 2 .^ str2double(z.log2copyratio);

S.READN.reads = 'not_available';
S.READT.reads = 'not_available';
