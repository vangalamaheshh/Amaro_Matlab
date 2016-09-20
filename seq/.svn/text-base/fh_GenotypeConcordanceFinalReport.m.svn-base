function fh_GenotypeConcordanceFinalReport(individual_name,tumor_metrics,normal_metrics,threshold,lane_blacklist)
% fh_GenotypeConcordanceFinalReport(individual_name,tumor_metrics,normal_metrics,threshold,lane_blacklist)
%
% Final step in the GenotypeConcordanceForIndividual pipeline
%   -- writes report.html
%   -- writes num.mixups.txt (number of lanes below threshold)
%
% Mike Lawrence 2010-01-21

if ~exist('threshold','var'), threshold = 95; end
if ~isnumeric(threshold), threshold = str2double(threshold); end

if ~exist('lane_blacklist','var'), lane_blacklist = 'none'; end
if strcmpi(lane_blacklist,'none')
  BL = [];
else
  BL = load_lines(lane_blacklist);
end

H = [...
  '<h3>Genotype Concordance Metrics for ' individual_name '</h3>'...
];

tn = {'tumor';'normal'};
fn = {tumor_metrics;normal_metrics};
no_data = false(2,1);

% LOAD RESULTS

L={};
for tni=1:2
   L{tni} = load_concordance_file(fn{tni});  
   if isempty(L{tni})
     H = [H '<p><font size=+1>No data for ' tn{tni} '</font>'];
     no_data(tni) = true;
   else 
     L{tni}.tn = repmat(tn(tni),slength(L{tni}),1);
   end
end
L = L(~no_data);

if ~isempty(L)
  L = concat_structs(L);

  L.is_blacklisted = is_blacklisted(L.lane,BL);    % apply blacklist

  L.pct_concordance = 100 * L.pct_concordance;  % make it ACTUALLY percent!
  
  tests = { 'HOMOZYGOUS_NON_REFERENCE','HOMOZYGOUS_REFERENCE','HETEROZYGOUS'};
  
  bad = (strcmp(L.category,tests{1}) & L.pct_concordance<threshold);
  mixups = find(bad & ~L.is_blacklisted);
  blacklisted_mixups = find(bad & L.is_blacklisted);
  
  % WRITE HTML REPORT
  
  if ~isempty(mixups)
    H = [H '<p><h2><font color="red">POSSIBLE MIXUPS: ' num2str(length(mixups)) '</font></h2>'];
  end
  if ~isempty(blacklisted_mixups)
    H = [H '<p><font color="blue">BLACKLISTED MIXUPS: ' num2str(length(blacklisted_mixups)) '</font>'];
  end
  if isempty(mixups) & isempty(blacklisted_mixups)
    H = [H '<p><h2>NO MIXUPS DETECTED</h2>'];
  end
  
  % lane tables

  L.is_blacklisted = nansub({'','+'},L.is_blacklisted+1);
  
  for t=1:3
    H = [H '<p><b>' tests{t} '</b>'];
    H = [H '<p><table class="its"><tr>'];
    
    flds1 = {'tn','lane','reference','non_reference','pct_concordance','is_blacklisted'};
    flds2 = {'T/N','lane','ref','nonref','%concordance','BL'};
    
    widths = [100 100 100 100 120 40];
    decplc = [nan nan 0   0   1   nan]; 
    commas = [nan nan 1   1   0   nan];
    rtjust = [0   0   1   1   1   1];
    
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
      if ~strcmp(L.category{i},tests{t}), continue; end
      H = [H '<tr>'];
      for f=1:length(flds1)
        H = [H '<td'];
        if rtjust(f)==1, H = [H ' align=right']; end
        H = [H '>'];
        if ismember(i,mixups), H = [H '<font color="#FF0000"><b>']; end   
        if ismember(i,blacklisted_mixups), H = [H '<font color="#0000FF"><b>']; end
        H = [H X{i,f}];
        if ismember(i,mixups) || ismember(i,blacklisted_mixups), H = [H '</font></b>']; end
        H = [H '</td>'];
      end
      H = [H '</tr>'];
    end
    
    H = [H '</table>'];
  end
end

H = [H '<p><hr>'...
  '<br><b>BL</b> = blacklisted lanes indicated with "+"'...
  '</html>'...
];

save_textfile(H,'report.html');

save_textfile(sprintf('%d\n',length(mixups)),'num.mixups.txt');
