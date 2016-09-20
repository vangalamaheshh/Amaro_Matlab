function H = parse_blast_output_file(filename)

% remove header lines and load file
tmpfile = [filename '_tmp'];
system(['cat ' filename ' | grep -v -P "^#" > ' tmpfile]);
f = fopen(tmpfile);
t = textscan(f,'%s%s%f%f%f%f%f%f%f%f%f%f');
fclose(f);
system(['rm ' tmpfile]);

H=[];
H.query = t{1};
H.hit = t{2};
H.pctid = t{3};
H.matchlen = t{4};
H.n_mm = t{5};
H.n_gaps = t{6};
H.qstart = t{7};
H.qend = t{8};
H.qstrand = repmat({'+'},slength(H),1);
H.hstart = t{9};
H.hend = t{10};
H.hstrand = H.qstrand;
H.E = t{11};
H.score = t{12};

% replace start>end with explicit strand indicators

idx = find(H.qstart>H.qend);
if ~isempty(idx)
  tmp = H.qstart(idx); H.qstart(idx) = H.qend(idx); H.qend(idx) = tmp;
  H.qstrand(idx) = repmat({'-'},length(idx),1);
end

idx = find(H.hstart>H.hend);
if ~isempty(idx)
  tmp = H.hstart(idx); H.hstart(idx) = H.hend(idx); H.hend(idx) = tmp;
  H.hstrand(idx) = repmat({'-'},length(idx),1);
end

