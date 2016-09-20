function H = parse_blast_result(result)

% if error/warning, blank out the results
if strncmp(result,'[blastall]',10), result=''; end

lines = split(result,char(10));   % split at newline
lines(strcmp('',lines)) = [];     % remove blank lines
lines(strncmp('#',lines,1)) = []; % remove comment lines

nlines = length(lines);

%lines,keyboard

H.query = cell(nlines,1);
H.hit = cell(nlines,1);
H.pctid = nan(nlines,1);
H.matchlen = nan(nlines,1);
H.n_mm = nan(nlines,1);
H.n_gaps = nan(nlines,1);
H.qstart = nan(nlines,1);
H.qend = nan(nlines,1);
H.qstrand = cell(nlines,1);
H.hstart = nan(nlines,1);
H.hend = nan(nlines,1);
H.hstrand = cell(nlines,1);
H.E = nan(nlines,1);
H.score = nan(nlines,1);
for i=1:nlines
  fields = split(lines{i},char(9));   % split at tabs
  H.query{i} = fields{1};
  H.hit{i} = fields{2};
  H.pctid(i) = str2double(fields{3});
  H.matchlen(i) = str2double(fields{4});
  H.n_mm(i) = str2double(fields{5});
  H.n_gaps(i) = str2double(fields{6});
  H.qstart(i) = str2double(fields{7});
  H.qend(i) = str2double(fields{8});
  H.hstart(i) = str2double(fields{9});
  H.hend(i) = str2double(fields{10});
  H.E(i) = str2double(fields{11});
  H.score(i) = str2double(fields{12});

  % replace start>end with explicit strand indicators

  if H.qstart(i)>H.qend(i)
    H.qstrand{i} = '-';
    [H.qstart(i),H.qend(i)] = exchange_vars(H.qstart(i),H.qend(i));
  else
    H.qstrand{i} = '+';
  end

  if H.hstart(i)>H.hend(i)
    H.hstrand{i} = '-';
    [H.hstart(i),H.hend(i)] = exchange_vars(H.hstart(i),H.hend(i));
  else
    H.hstrand{i} = '+';
  end

end


end
