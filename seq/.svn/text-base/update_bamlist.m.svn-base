function B = update_bamlist(B)
% B = update_bamlist(B)
%
% given B, with required field B.bam, a list of BAM filenames,
% does the following:
%
% (1) verifies that B.bam exists, adding B.ok
% (2) if it doesn't exist, tries to rescue by incrementing the version number
% (3) if it does exist, looks up the "analysis_files.txt" to find out:
%     -- date, datenum (earliest and latest lanes)
%     -- build
%     -- baitset
%     -- aligner

require_fields(B,'bam');
nb = slength(B);

B.ok = false(nb,1);
z = repmat({'-'},nb,1);
B.build = z;
B.aligner = z;
B.date_earliest = z;
B.datenum_earliest = nan(nb,1);
B.date_latest = z;
B.datenum_latest = nan(nb,1);

% fields to extract from analysis_files.txt
flds = {'bait_set','analysis_type','initiative'};
flds2 = {'baitset','analysis_type','initiative'};
for i=1:length(flds2)
  B.(flds2{i}) = z;
end

% SURVEY THE BAMS
for i=1:nb
  if ~mod(i,10), fprintf('%d/%d ',i,nb); end
  found = false;
  while(true)
    b = parse(B.bam{i},'^(.*)+/([^/]+)$',{'path','bam'});
    b = parse_in(b,'path','^(.*)+/v(\d+)$',{'basepath','ver'},2);
    if exist(B.bam{i})
      found = true;
      break
    else
      if isempty(b.basepath) || ~exist(b.basepath{1},'dir'), break; end   %% doesn't appear to have a version structure
      d=[]; d.bam = direc([b.basepath{1} '/v*/' b.bam{1}]);
      if isempty(d.bam), break; end   %% no versions available at all
      d = parse_in(d,'bam','^(.*)+/([^/]+)$',{'path','bamname'});
      d = parse_in(d,'path','^(.*)+/v(\d+)$',{'basepath','ver'},2);
      d = sort_struct(d,'ver',-1);
      B.bam{i} = d.bam{1};
      continue
    end
  end
  B.ok(i) = found;
  if ~found, continue; end  %% could not find this BAM
  % process BAM header
  try
    bam = open_BAM_file(B.bam{i});
    h = char(bam.getFileHeader.getTextHeader.toString);
    used_maq = contains(h,'MAQ') || contains(h,'maq');
    used_bwa = contains(h,'BWA') || contains(h,'bwa');
    used_hg18 = contains(h,'hg18') || contains(h,'Homo_sapiens_assembly18');
    used_hg19 = contains(h,'hg19') || contains(h,'Homo_sapiens_assembly19');
    if used_maq && used_bwa
      B.aligner{i} = 'maq+BWA';
    elseif used_maq
      B.aligner{i} = 'maq';
    elseif used_bwa
      B.aligner{i} = 'BWA';
    else
      B.aligner{i} = '???';
    end
    if used_hg18 && used_hg19
      B.build{i} = 'hg18+hg19';
    elseif used_hg18
      B.build{i} = 'hg18';
    elseif used_hg19
      B.build{i} = 'hg19';
    else
      B.build{i} = '???';
    end
    d = parse(grep('^@RG',text_to_lines(h)),'DT:(\d\d\d\d)-(\d\d)-(\d\d)',{'yr','mo','day'});
    d.date = stringsplice([d.yr d.mo d.day],'-');
    d.datenum = datenum(d.date);
    d = sort_struct(d,'datenum');
    B.date_earliest{i} = d.date{1};
    B.datenum_earliest(i) = d.datenum(1);
    B.date_latest{i} = d.date{end};
    B.datenum_latest(i) = d.datenum(end);
    bam.close();
  catch me
  end
  % process analysis_files.txt
  try
    af = [b.path{1} '/analysis_files.txt'];
    if exist(af,'file')
      params=[]; params.lowercase_fieldnames = true;
      A = load_struct(af,params);
      if ~isfield(A,'bait_set') && isfield(A,'bait_intervals')
        A.bait_set = regexprep(A.bait_intervals,'.*HybSelOligos/([^/]+)/.*$','$1');
      end
      if ~isfield(A,'analysis_type') && isfield(A,'library_type')
        A = rename_field(A,'library_type','analysis_type');
      end
      % extract desired fields
      for f=1:length(flds)
        if isfield(A,flds{f})
          tmp = getfield(A,flds{f});
          tmp = tmp(~strcmp('',tmp));
          B.(flds2{f}){i} = concat(unique(tmp),',');
        end
      end
    end
  catch me
  end
end
fprintf('Done.\n');






  
