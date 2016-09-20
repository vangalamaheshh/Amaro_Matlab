function H = dRanger_write_report(individual,circospng,L,HIGH_CONFIDENCE_THRESHOLD)
% dRanger_write_report(individual,circospng,L)
%
% Mike Lawrence 2010

if ~exist('HIGH_CONFIDENCE_THRESHOLD','var')
  HIGH_CONFIDENCE_THRESHOLD= 4;
end

H = ['<h3>dRanger results for ' individual '</h3>'];

inter = ~strcmp(L.chr1,L.chr2);
span = L.span;
if ~isnumeric(span), span = str2double(span); end
small = (span<10000);
med = (span>=10000 & span<1e6);
large = (span>=1e6);
hc = L.score>=HIGH_CONFIDENCE_THRESHOLD;

circospng = regexprep(circospng,'.*/([^/])+','$1');

H = [H '<p><table><tr><td width=400>'...
];

if any(~hc)
  H = [H ...
  '<b>High-confidence somatic rearrangements<br>(score &ge;' num2str(HIGH_CONFIDENCE_THRESHOLD) ' <br>shown on CIRCOS plot)</b>'...
  ];
else
  H = [H ...
  '<b>Somatic rearrangements<br>(score &ge;' num2str(HIGH_CONFIDENCE_THRESHOLD) ' <br>shown on CIRCOS plot)</b>'...
  ];
end

H = [H ...
  '<table><tr><td width=250 align=right>Interchromosomal:</td><td width=40 align=right>' num2str(sum(inter&hc)) '</td></tr>'...
  '<tr><td align=right>Long-range (1Mb+):</td><td align=right>' num2str(sum(~inter&large&hc)) '</td></tr>'...
  '<tr><td align=right>Mid-range (10Kb&ndash;1Mb):</td><td align=right>' num2str(sum(~inter&med&hc)) '</td></tr>'...
  '<tr><td align=right>Local (<10Kb):</td><td align=right>' num2str(sum(~inter&small&hc)) '</td></tr>'...
  '<tr><td align=right>Total:</td><td align=right>' num2str(sum(hc)) '</td></tr></table>'...
];

if any(~hc)
  H = [H ...
       '<p><b>All candidate somatic rearrangements</b>'...
       '<table><tr><td width=250 align=right>Interchromosomal:</td><td width=40 align=right>' num2str(sum(inter)) '</td></tr>'...
       '<tr><td align=right>Long-range (1Mb+):</td><td align=right>' num2str(sum(~inter&large)) '</td></tr>'...
       '<tr><td align=right>Mid-range (10Kb&ndash;1Mb):</td><td align=right>' num2str(sum(~inter&med)) '</td></tr>'...
       '<tr><td align=right>Local (<10Kb):</td><td align=right>' num2str(sum(~inter&small)) '</td></tr>'...
       '<tr><td align=right>Total:</td><td align=right>' num2str(slength(L)) '</td></tr></table>'...
      ];
end

H = [H ...
   '</td>'...
   '<td width=600>'...
   '<img src="' circospng '" width=500>'...
   '</td></tr></table>' ...
];

H = [H '<hr>'];

% lane table
H = [H '<table class="its"><tr>'];

flds1 = {'num','chr1','str1','pos1','chr2','str2','pos2','tumreads','normreads',...
  'class','span',...
  'site1','site2','fusion','quality','score','BPresult','BPsomaticratio','T_lenhomology','T_lenforeign'};
flds2 = {'#','chr1','str1','pos1','chr2','str2','pos2','#T','#N',...
  'class','span',...
  'site1','site2','fusion details','qual','score','BP','BPscore','homology length','nontempleted sequence length'};

widths = [40 30  25  120 30  25  120   30 30   75   75  170 170 140 40 40 25 50 50 50];
decplc = [ 0 nan nan  0  nan nan  0    0  0   nan   0   nan nan nan  2  1 0 1 1 1];
commas = [ 0 nan nan  1  nan nan  1    0  0   nan   1   nan nan nan  0  0 0 0 0 0];
justfy = [-1 -1 -1   -1  -1 -1   -1   -1 -1   -1   -1   -1  -1  nan -1 -1 -1 -1 -1 -1];

L.str1 = regexprep(L.str1,'\-','&ndash;');
L.str2 = regexprep(L.str2,'\-','&ndash;');

if ~isnumeric(L.BPresult), L.BPresult = str2double(L.BPresult); end
L.BPresult = nansub({'germ';' ';'+'},L.BPresult+2,'&mdash;');

L.fusion_priority = zeros(slength(L),1);
idx = grep('in frame|Protein fusion: mid-exon',L.fusion,1);
L.fusion_priority(idx) = 1;

X = cell(slength(L),length(flds1));
for f=1:length(flds1)
  tmp = getfield(L,flds1{f});
  if isnumeric(tmp)
    for i=1:length(tmp)
      if isnan(tmp(i))
        X{i,f} = '&mdash;';
      else
        X{i,f} = sprintf(['%0.' num2str(decplc(f)) 'f'],tmp(i));
        if commas(f)==1, X{i,f} = comma_format(X{i,f}); end
      end
    end
  else
    X(:,f) = tmp;
  end
end

for f=1:length(flds1)
  H = [H '<th'];
  if justfy(f)==1, H = [H ' align=right']; end
  if justfy(f)==-1, H = [H ' align=left']; end
  H = [H ' width="' num2str(widths(f)) '">' flds2{f} '</th>'];
end
H = [H '</tr>'];

fidx1 = find(strcmp(flds1,'fusion'));
fidx2 = [find(strcmp(flds1,'pos1'));find(strcmp(flds1,'pos2'))];
for i=1:slength(L)
  H = [H '<tr>'];
  for f=1:length(flds1)
    H = [H '<td'];
    if justfy(f)==1, H = [H ' align=right']; end
    if justfy(f)==-1, H = [H ' align=left']; end
    H = [H '>'];
    redflag = (f==fidx1 && L.fusion_priority(i)>=1);
    approxflag = (ismember(f,fidx2) && L.approxflag(i));
    if redflag, H = [H '<font color="red"><b>']; end     
    if approxflag, H = [H '<font color="gray">[~ ']; end
    H = [H X{i,f}];
    if redflag, H = [H '</b></font>']; end
    if approxflag, H = [H ']</font>']; end
    H = [H '</td>'];
  end
  H = [H '</tr>'];
end

H = [H '</table><p><hr>'...
  '<br><b>chr1 str1 pos1 site1</b> = genomic position of one side of the rearrangement'...
  '<br><b>chr2 str2 pos2 site2</b> = genomic position of the other side of the rearrangement'...
  '<p><u>Note:</u>  <b>pos1</b> and <b>pos2</b> are BreakPointer-refined coordinates except when shown in brackets.'...
  '<p><u>Note:</u>  <b>site1</b> and <b>site2</b> are listed here in ascending genomic order; there is no physical significance'...
  '<br>to this ordering.  The two sites represent the two 5'' ends of the chimeric DNA molecules'...
  '<br>that were sequenced and found to align to nonadjacent genomic loci.  Note that inversion events'...
  '<br>are characterized by having <b>str1</b> = <b>str2</b>.  Non-inversion events have <b>str1</b> &ne; <b>str2</b>.'...
  '<p><b>#T</b> = number of supporting readpairs in the tumor'...
  '<br><b>#N</b> = number of supporting readpairs in the normal (typically zero for somatic events)'...
  '<br><b>class</b> = description based on simplistic interpretation of relative positions and strands'...
  '<br><b>span</b> = for intrachromosomal events, the distance between Site1 and Site2'...
  '<br><b>fusion details</b> = gives details about rearrangements where both endpoints are in transcribed regions'...
  '<br><b>qual</b> = quality of candidate rearrangement from 1 (high-quality) to 0 (low-quality, likely artifactual)'...
  '<br><b>score</b> = computed by multiplying <b>#T</b> &times; <b>qual</b>'...
  '<br><b>BP</b> = BreakPointer result: &mdash; = not tried; blank = inconclusive; germ = germline; + = confirmed somatic'...
  '</html>'...
];

save_textfile(H,'report.html');
