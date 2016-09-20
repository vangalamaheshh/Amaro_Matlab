function convert_target_list(infile,outfile)
% convert_target_list(infile,outfile)
%
% Converts a list of hybrid-selection capture targets to the format
% used for the CancerGenomeAnalysis sequencing analysis pipeline.
% (removes header, imputes gene names, and adds GC content).
%
% Mike Lawrence 2009-06-09
%
% infile should have lines like the following (tab-delimited):
% chr1    939225  939723  +       target_67
% chr1    945414  945618  +       target_68
% chr1    947442  947707  +       target_69
% chr1    960518  960569  +       target_70
% chr1    965906  966125  +       target_71
% ...
% chrY    21150880        21150967        +       target_185720
% chrY    21153862        21153969        +       target_185721
% chrY    21155746        21155800        +       target_185722

% only the first three columns are used.
% lines not beginning with "chr" are ignored.
%
% outfile will have the following format (tab-delimited):
% the columns (no header line) are: <gene> <chr(1-24)> <start> <end> <gc_content>
% TNFRSF4 1       1136798 1136868 6.619718e-01
% TNFRSF4 1       1136947 1137075 7.441860e-01
% TNFRSF4 1       1137185 1137381 7.055838e-01
% TNFRSF4 1       1137881 1137947 6.716418e-01
% TNFRSF4 1       1138235 1138336 6.666667e-01

build = 'hg18'

try

fprintf('Reading %s\n',infile);
L = load_lines(infile); origlen = length(L);
idx = grep('^chr',L,1);
fprintf('  File contains %d lines; using the %d lines that begin with "chr".\n',origlen,length(idx));
skipped = setdiff(min(idx):max(idx),idx);
if ~isempty(skipped)
  fprintf('  WARNING:  File contains %d internal lines that had to be skipped:\n',length(skipped));
  disp(L(skipped));
end
L = L(idx);
try T = parse(L,'^(chr[^\t]*)\t(\d*)\t(\d*)',{'chrom','start','end'});
catch me, disp(me.message), error('Problem parsing lines!'); end
T.chr = convert_chr(T.chrom);
T.start = str2double(T.start); T.end = str2double(T.end);
if any(isnan(T.chr))
  fprintf(['WARNING:  Targets outside chr1...22XY (e.g. "chrM", "chrN_random") are not currently handled properly.\n'...
  'Infile contains %d of these targets; they will be omitted from outfile:\n'],sum(isnan(T.chr)));
  disp(L(isnan(T.chr)));
  T = reorder_struct(T,~isnan(T.chr));
end
nt = slength(T);

R = load_refseq(build);
R.chr = convert_chr(R.chr);

fprintf('Imputing gene names\n');
T.gene = cell(nt,1);

fprintf('  pass 1/4: if hits only one gene''s exon(s)\n');
for i=1:nt
  if ~mod(i,10000), fprintf('%d/%d ',i,nt); end
  idx = find(R.chr==T.chr(i) & R.tx_start <= T.end(i) & R.tx_end >= T.start(i));
  g = unique(R.gene(idx));
  if isempty(g)  fprintf_silent('%d. None found\n',i);
  elseif length(g)>1
    exonic = []; for j=1:length(idx), for e=1:R.n_exons(idx(j))
        if R.exon_starts{idx(j)}(e) <= T.end(i) & R.exon_ends{idx(j)}(e) >= T.start(i)
           exonic = [exonic; idx(j)];
    end,end,end
    g_exonic = unique(R.gene(exonic));
    if isempty(g_exonic)  fprintf_silent('%d. None exonic found\n',i);
    elseif length(g_exonic)>1, fprintf_silent('%d. Multiple exonic found:',i);
      for j=1:length(g_exonic), fprintf_silent(' %s',g_exonic{j}); end, fprintf_silent('\n');
    else T.gene{i}=g_exonic{1}; end
  else T.gene{i}=g{1}; end
end

fprintf('  %d/%d targets assigned to genes\n',sum(~cellfun('isempty',T.gene)),nt);
fprintf('\n  pass 2/4: if hits exons of adjacent gene on list\n');
while(1)
n = sum(cellfun('isempty',T.gene));
for i=1:nt
  if ~isempty(T.gene{i}), continue; end
  idx = find(R.chr==T.chr(i) & R.tx_start <= T.end(i) & R.tx_end >= T.start(i));
  g = unique(R.gene(idx));
  if length(g)>1
    adjacent = {};
    if i>1, adjacent = [adjacent; T.gene{i-1}]; end
    if i<nt, adjacent = [adjacent; T.gene{i+1}]; end
    g_adjacent = intersect(g,adjacent);
    if isempty(g_adjacent)  fprintf_silent('%d. None adjacent found\n',i);
    elseif length(g_adjacent)>1, fprintf_silent('%d. Multiple adjacent found:',i);
      for j=1:length(g_adjacent), fprintf_silent(' %s',g_adjacent{j}); end, fprintf_silent('\n');
    else T.gene{i}=g_adjacent{1}; end
  end
end
n_new = sum(cellfun('isempty',T.gene));
if (n_new==n) break; end
end

fprintf('  %d/%d targets assigned to genes\n',sum(~cellfun('isempty',T.gene)),nt);
fprintf('\n  pass 3/4: first overlapping gene\n');
for i=1:nt
  if ~isempty(T.gene{i}), continue; end
  idx = find(R.chr==T.chr(i) & R.tx_start <= T.end(i) & R.tx_end >= T.start(i));
  g = unique(R.gene(idx));
  if length(g)>1, T.gene{i}=g{1}; end
end

fprintf('  %d/%d targets assigned to genes\n',sum(~cellfun('isempty',T.gene)),nt);
fprintf('\n  pass 4/4: nearby genes\n');
for margin=[1000,2000,5000,10000,20000,50000,100000,200000,500000,10000000]
  for i=1:nt
    if ~isempty(T.gene{i}), continue; end
    idx = find(R.chr==T.chr(i) & R.tx_start <= T.end(i)+margin & R.tx_end >= T.start(i)-margin);
    g = unique(R.gene(idx));
    if length(g)>1, T.gene{i}=g{1}; end
  end
end
fprintf('  %d/%d targets assigned to genes\n',sum(~cellfun('isempty',T.gene)),nt);

nempty = sum(cellfun('isempty',T.gene));
if nempty>0
  fprintf('Unable to assign %d targets to genes--> marking these "Unknown"\n',nempty);
  T.gene = fillblanks(T.gene,'Unknown');
end

% GC content
fprintf('\nAnnotating GC content\n');
T.gc = nan(nt,1);
for i=1:nt
  if ~mod(i,10000), fprintf('%d/%d ',i,nt); end
  ref = upper(genome_region(T.chr(i),T.start(i),T.end(i),build));
  gc = nan(length(ref),1);
  gc(ref=='A' | ref=='T') = 0;
  gc(ref=='C' | ref=='G') = 1;
  T.gc(i) = nanmean(gc);
end

% save file
fprintf('\nSaving file %s\n',outfile);
T = keep_fields(T,{'gene','chr','start','end','gc'});
save_struct(T,outfile,'no_headers');
fprintf('Done!\n');

catch me, excuse(me); end
