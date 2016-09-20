function H = convert_CNQC_report_to_html(L,Z,individual_name)
% Mike Lawrence 2009-10-29

H = [...
  '<html><head><title>CopyNumberQC for ' individual_name '</title></head>' ...
  '<h2>CopyNumberQC for ' individual_name '</h2>'
];

mixups = find(L.is_mixup & ~L.is_blacklisted);
blacklisted_mixups = find(L.is_mixup & L.is_blacklisted);
if ~isempty(mixups)
  H = [H '<p><h2><font color="red">POSSIBLE MIXUPS: ' num2str(length(mixups)) '</font></h2>'];
end
if ~isempty(blacklisted_mixups)
  H = [H '<p><font color="blue">BLACKLISTED MIXUPS: ' num2str(length(blacklisted_mixups)) '</font>'];
end
if isempty(mixups) & isempty(blacklisted_mixups)
  H = [H '<p><h2>NO MIXUPS DETECTED</h2>'];
end

% plot

H = [H '<hr><h2>Copy number by lane</h2>'];

if ~any(L.enoughreads(strcmp('tumor',L.tn)))
  H = [H '<p>No tumor lanes had enough reads to analyze.'];
end

if ~any(L.enoughreads(strcmp('normal',L.tn)))
  H = [H '<p>No normal lanes had enough reads to analyze.'];
end

H = [H '<p><img src="' individual_name '_CopyNumberQC.png">'];

% plot legend

%if ~isempty(mixups), H = [H '<br><font size=+2><b>*</b></font> = lane flagged as mixup']; end
%if any(L.is_blacklisted) H = [H '<br><font color="000000" size=+2><b>*</b></font> = blacklisted lane']; end

if ~isempty(grep('tumor',L.tn(mixups))), H = [H '<br><font color="FF0000" size=+2><b>*</b></font> = tumor-lane mixup']; end
if ~isempty(grep('normal',L.tn(mixups))), H = [H '<br><font color="FF8888" size=+2><b>*</b></font> = normal-lane mixup'];end
if ~isempty(grep('tumor',L.tn(blacklisted_mixups))), H = [H '<br><font color="0000FF" size=+2><b>*</b></font> = blacklisted tumor lane']; end
if ~isempty(grep('normal',L.tn(blacklisted_mixups))), H = [H '<br><font color="8888FF" size=+2><b>*</b></font> = blacklisted normal lane'];end

if any(~isnan(Z(:,1))), H = [H  '<br><b>S<sub>T</sub></b> = array-based segfile for tumor']; end
if any(~isnan(Z(:,2))), H = [H  '<br><b>S<sub>N</sub></b> = array-based segfile for normal']; end

% lane table

L.is_blacklisted = nansub({'';'+'},L.is_blacklisted+1);

H = [H '<hr><h2>Lane details</h2>'];

if any(~L.enoughreads)
  H = [H '<p>Low-counts lanes are greyed out and do not appear in the above plot.'];
end

H = [H '<p><table class="its"><tr>'];

flds1 = {'tn','readgroup','flowcell_lane','library','baitset','nreads',...
  'noise','normalness','tumorseg_corr','tumormed_corr','judgement','is_blacklisted'};
flds2 = {'T/N','rgrp','flowcell_lane','library','baitset','#reads',...
  'noise','N-ness','Tcorr1','Tcorr2','judgement','BL'};

widths = [68  80  200 120 20  85    50 50 50 50 210 20];
decplc = [nan nan nan nan nan 0     2  2  2  2  nan nan];
commas = [nan nan nan nan nan 1     0  0  0  0  nan nan];
rtjust = [0   0   0   0   0   1     1  1  1  1  1   1];

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
  H = [H '<tr>'];
  for f=1:length(flds1)
    H = [H '<td'];
    if rtjust(f)==1, H = [H ' align=right']; end
    H = [H '>'];
    if ~L.enoughreads(i)
      H = [H '<font color=grey>'];
    end
    if ismember(i,mixups)
      if strcmp(L.tn{i},'tumor')
        H = [H '<font color="#FF0000"><b>'];
      else
        H = [H '<font color="#FF8888"><b>'];
      end
    end
    if ismember(i,blacklisted_mixups)
      if strcmp(L.tn{i},'tumor')
        H = [H '<font color="#0000FF"><b>'];
      else
        H = [H '<font color="#8888FF"><b>'];
      end
    end
    H = [H X{i,f}];
    if ismember(i,mixups) || ismember(i,blacklisted_mixups), H = [H '</b></font>']; end
    if ~L.enoughreads(i), H = [H '</font>']; end
    H = [H '</td>'];
  end
  H = [H '</tr>'];
end

% table legend

H = [H '</table><p><hr>'...
  '<br><b>noise</b>  = noise (>0.40 = ''extremely noisy'')'...
  '<br><b>N-ness</b> = normalness (&lt;0.8 = tumor, &gt;0.9=normal)'...
  '<br><b>Tcorr1</b> = correlation to tumor array (if available)'...
  '<br><b>Tcorr2</b> = correlation to median of tumor lanes'...
  '<br><b>BL</b> = blacklisted lanes indicated with "+"'...
  '</html>'...
];
