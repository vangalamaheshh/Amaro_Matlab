function fh_global_coverage_stats_report(individual_name, filterdir, gczdir, outname)
% fh_global_coverage_stats_report(individual_name, filterdir, gczdir)
%
% individual_name = e.g. GBM-0188
%
% filterdir = should be the full path of a directory containing the 25 filtering files:
%                     chrN.txt for N=1-24
%                            with one line per basepair with an integer indicating the filtering category
%                     categs.txt listing the names of the filtering categories (with header line)
%             e.g. /xchip/cga/home/lawrence/reference/genome_filters/zone
%             or else 'none' for no category-based filtering
%
% gczdir = output directory of GetCoverageStatsUsingFilter, containing:
%                     chrN.txt for N=1-24
%
% outputs:
%           <outname> = computer-readable
%           report.html = human-readable
%
% Mike Lawrence 2009

if nargin~=4
  error(['input parameters required:\n'...
        '   individual_name, filterdir, gczdir, outname\n']);
end

% load category names
F = load_struct([filterdir '/categs.txt'],'%f%s');

% integrate statistics across chromosomes

S = cell(24,1);
for c=1:24
  statfile = [gczdir '/chr' num2str(c) '.txt'];
  S{c} = load_struct(statfile,'%f%f%f%f%f%f',0);
end
S = rename_fields(concat_structs(S),{'col1','col2','col3','col4','col5','col6'},...
                  {'chr','categ','terr','tseqbp','nseqbp','callable'});
[categs ui uj] = unique(S.categ);
ncateg = length(categs);
L = [];
for j=1:ncateg
  L.categ_num{j,1} = num2str(categs(j));
  L.categ_name{j,1} = F.name{F.num==categs(j)};
  L.terr(j,1) = sum(S.terr(uj==j));
  L.tseqbp(j,1) = sum(S.tseqbp(uj==j));
  L.nseqbp(j,1) = sum(S.nseqbp(uj==j));
  L.tdepth(j,1) = L.tseqbp(j) / L.terr(j);
  L.ndepth(j,1) = L.nseqbp(j) / L.terr(j);
  L.callablebp(j,1) = sum(S.callable(uj==j));
  L.callablepct(j,1) = 100 * L.callablebp(j) / L.terr(j);
end
% add row of summary statistics
j = ncateg+1;
L.categ_num{j,1} = '-';
L.categ_name{j,1} = 'Total';
L.terr(j,1) = sum(L.terr(1:j-1));
L.tseqbp(j,1) = sum(L.tseqbp(1:j-1));
L.nseqbp(j,1) = sum(L.nseqbp(1:j-1));
L.tdepth(j,1) = L.tseqbp(j) / L.terr(j);
L.ndepth(j,1) = L.nseqbp(j) / L.terr(j);
L.callablebp(j,1) = sum(L.callablebp(1:j-1));
L.callablepct(j,1) = 100 * L.callablebp(j) / L.terr(j);

if isempty(L)
  error('WEIRD PROBLEM:  Empty L during fh_global_coverage_stats_report for %s!\n',individual_name);
end

% write textfile output

save_struct(L,outname);

% write html report

L.categ_num{end} = '';  % don't display first-column "-" in html report

H = [...
  '<h3>Global Coverage Statistics for ' individual_name '</h3>'...
  '<hr>'...
];

% lane table
H = [H '<table class="its"><tr>'];

flds1 = {'categ_num','categ_name','terr','tseqbp','nseqbp',...
  'tdepth','ndepth','callablebp','callablepct'};
flds2 = {'#','category','territory','sequenced bp<br>(Tumor)','sequenced bp<br>(Normal)'...
  'depth<br>(Tumor)','depth<br>(Normal)','callable bp','%callable'};

widths = [30  80  130 140 140 70  70 130 90];
decplc = [nan nan 0   0   0   1   1  0   1];
commas = [nan nan 1   1   1   0   0  1   0];
rtjust = [1   1   1   1   1   1   1  1   1];

X = cell(slength(L),length(flds1));
for f=1:length(flds1)
  tmp = getfield(L,flds1{f});
  if isnumeric(tmp)
    for i=1:length(tmp)
      X{i,f} = sprintf(['%0.' num2str(decplc(f)) 'f'],tmp(i));
      if commas(f)==1, X{i,f} = comma_format(X{i,f}); end
    end
  else
    X(:,f) = tmp;
  end
end

for f=1:length(flds1)
  H = [H '<th'];
  if rtjust(f)==1, H = [H ' align=right']; end
  H = [H ' width="' num2str(widths(f)) '">' flds2{f} '</th>'];
end
H = [H '</tr>'];

for i=1:slength(L)
  % break for last row
  if i==slength(L), H = [H '<tr></tr>']; end
  H = [H '<tr>'];
  for f=1:length(flds1)
    H = [H '<td'];
    if rtjust(f)==1, H = [H ' align=right']; end
    H = [H '>'];
    if i==slength(L) && f<3, H = [H '<b>']; end
    H = [H X{i,f}];
    if i==slength(L) && f<3, H = [H '</b>']; end
    H = [H '</td>'];
  end
  H = [H '</tr>'];
end

H = [H '</table><p><hr>'...
  '<br><b>territory</b> = total genomic territory in this category'...
  '<br><b>sequenced bp</b> = total number of sequenced basepairs in this category'...
  '<br><b>depth</b> = average sequencing depth achieved in this category'...
  '<br><b>callable bp</b> = total number of callable basepairs (14/8) in this category'...
  '<br><b>%callable</b> = percent of territory callable in this category'...
  '</html>'...
];

save_textfile(H,'report.html');
