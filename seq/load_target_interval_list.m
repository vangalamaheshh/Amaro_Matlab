function T = load_target_interval_list(infile)
% Mike Lawrence 2010

x = load_lines(infile);
% discard metadata lines and comment lines
x = grep('^[^#@]',x);

% parse into chr,start,end
T = parse(x,'^(\S*)\t(\d*)\t(\d*)\t.*',{'chrname','start','end'});
T.chr = convert_chr(T.chrname);
T = make_numeric(T,{'start','end'});
