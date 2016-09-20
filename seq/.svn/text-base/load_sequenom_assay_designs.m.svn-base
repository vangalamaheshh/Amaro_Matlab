function A = load_sequenom_assay_designs(fname)
demand_file(fname);
A.line = load_lines(fname);
A = parse_in(A,'line','^(.*)[\t,].*$',{'name'});
A = parse_in(A,'name','^(chr[\dXY]+)_(\d+)(_gI_|_gD_|_)(\d+)$',{'chr','assay_pos','aux','num'},[2 4]);
A.aux = regexprep(A.aux,'_','');
A.chr = convert_chr(A.chr);
A.pos = A.assay_pos;
idx = grep('D',A.aux,1); A.pos(idx) = A.pos(idx) + 1;    % adjust deletion coordinates to match our convention
