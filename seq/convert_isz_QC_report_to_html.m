function H = convert_isz_QC_report_to_html(R,individual_name)
% Mike Lawrence 2009-10-15

H = [...
  '<html><head><title>InsertSizeQC for ' individual_name '</title></head>' ...
  '<h1>InsertSizeQC for ' individual_name '</h1>'
];

mixups = find(R.is_mixup & ~R.is_blacklisted);
blacklisted_mixups = find(R.is_mixup & R.is_blacklisted);

if all(isnan(R.adjmean))
  H = [H '<p><h2><font color="red">NO DATA</font></h2>'];
else
  if ~isempty(mixups)
    H = [H '<p><h2><font color="red">POSSIBLE MIXUPS: ' num2str(length(mixups)) '</font></h2>'];
  end
  if ~isempty(blacklisted_mixups)
    H = [H '<p><font color="blue">BLACKLISTED MIXUPS: ' num2str(length(blacklisted_mixups)) '</font>'];
  end
  if isempty(mixups) & isempty(blacklisted_mixups)
    H = [H '<p><h2>NO MIXUPS DETECTED</h2>'];
  end
end



% insert plots
H = [H ...
  '<hr><p><h2>Insert size distributions</h2><br><img src="' individual_name '_insert_size_distributions.png">' ...
  '<hr><p><h2>Insert size by flowcell</h2><br><img src="' individual_name '_insert_size_by_flowcell.png">' ...
  '<hr><p><h2>Peak mean/width scatterplot</h2><br><img src="' individual_name '_insert_size_scatterplot.png">' ...
];

% lane table
H = [H '<hr><p><h2>Lane details</h2><br><table class="its"><tr>'];

R.is_blacklisted = nansub({'','+'},double(R.is_blacklisted)+1);

flds1 = {'tn','readgroup','flowcell_lane','library','nreads',...
  'adjmean','width','adjmean_pctdev','width_pctdev','judgement','is_blacklisted'};

flds2 = {'T/N','rgrp','flowcell_lane','library','#reads',...
  'adjmean','width','%&Delta;am','%&Delta;w','judgement','BL'};

widths = [70 50 180 120 90 70 60 50 50 120 20];
decplc = [nan nan nan nan 0 2 1 1 1 nan nan];

X = cell(slength(R),length(flds1));
for f=1:length(flds1)
  tmp = getfield(R,flds1{f});
  if isnumeric(tmp)
    for i=1:length(tmp), X{i,f} = sprintf(['%0.' num2str(decplc(f)) 'f'],tmp(i)); end
  else
    X(:,f) = tmp;
  end
end

for f=1:length(flds1)
  H = [H '<th width="' num2str(widths(f)) '">' flds2{f} '</th>'];
end
H = [H '</tr>'];

for i=1:slength(R)
  H = [H '<tr>'];
  for f=1:length(flds1)
    H = [H '<td>'];
    if ismember(i,mixups)
      if strcmp(R.tn{i},'tumor')
        H = [H '<font color="#FF0000"><b>'];
      else
        H = [H '<font color="#FF8888"><b>'];
      end
    end
    if ismember(i,blacklisted_mixups)
      if strcmp(R.tn{i},'tumor')
        H = [H '<font color="#0000FF"><b>'];
      else
        H = [H '<font color="#8888FF"><b>'];
      end
    end   
    H = [H X{i,f}];
    if ismember(i,mixups) || ismember(i,blacklisted_mixups), H = [H '</b></font>']; end
    H = [H '</td>'];
  end
  H = [H '</tr>'];
end

H = [H '</table><p><hr>'...
  '<br><b>adjmean</b> = adjusted mean insert size'...
  '<br><b>width</b> = peak width at half height'...
  '<br><b>%&Delta;am</b> = percent deviation from median adjmean'...
  '<br><b>%&Delta;w</b> = percent deviation from median width'...
  '<br><b>BL</b> = blacklisted lanes indicated with "+"'...
  '</html>'...
];

