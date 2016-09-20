function P = crawl_picard(rikerdir)
% Mike Lawrence 2009-09-15

if nargin==0, error('please specify rikerdir'); end

P = [];
P.project = [];
P.sample = [];
P.version = [];

flds = {...
      'initiative';...
      'work_request_id';...
      'sample_barcode';...
      'sample_alias';...
      'paired_run';...
      'library_name';...
      'flowcell_lane';...
      'bait_set';...
      'analysis_type';...
      'reference_sequence';...
      'target_intervals';...
      'bait_intervals';...
};
F = {};
for f=1:length(flds)
  P = setfield(P,flds{f},[]);
end

P.numlanes = [];
P.bampath = [];
P.bamsize = [];
P.bamtime = [];
P.bamtimenum = [];

P.finished = [];

idx=0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CRAWL

picdns = {'/seq/picard/by_project';'/seq/picard_aggregation'};
for q=1:length(picdns)
 picdn = picdns{q};
 picd = dir(picdn);
 for i=1:length(picd)
  proj = picd(i).name;
  if proj(1)=='.', continue; end
  projdn = [picdn '/' proj];
  projd = dir(projdn);
  for j=1:length(projd)
    samp = projd(j).name;
    if samp(1)=='.', continue; end
    sampdn = [projdn '/' samp];
    sampd = dir(sampdn);
    for k=1:length(sampd)
      ver = sampd(k).name;
      if ver(1)=='.', continue; end
      if strcmpi(ver,'current'), continue; end    % ignore "current" shortcut
      verdn = [sampdn '/' ver];
      idx=idx+1;
      P.project{idx,1} = proj; 
      P.sample{idx,1} = samp;
      P.version{idx,1} = ver;
      for f=1:length(flds), F{idx,f} = '?'; end
      analysisfn = [verdn '/analysis_files.txt'];
      processed_afile = false;
      if exist(analysisfn,'file')
        params=[]; params.lowercase_fieldnames = true;
        try
          A = load_struct(analysisfn,params);
          if ~isempty(A) && slength(A)>0
            % remove comment lines
            tmp = fieldnames(A);
            tmp = getfield(A,tmp{1});
            tmp = grep('^#.*',tmp,1);
            A = reorder_struct(A,setdiff(1:slength(A),tmp));
            % extract desired fields
            for f=1:length(flds)
              if isfield(A,flds{f})
                tmp = getfield(A,flds{f});
                tmp = tmp(~strcmp('',tmp));
                F{idx,f} = concat(unique(tmp),',');
              end
            end
            P.numlanes(idx,1) = slength(A);
            processed_afile = true;
          end         
        catch me
          % failed to process analysis_files.txt
        end
        if ~processed_afile
          for f=1:length(flds)
            F{idx,f} = '---';
          end
          P.numlanes(idx,1) = 0;
        end
      end
      bamn = [verdn '/' samp '.bam'];
      if ~exist(bamn,'file')
        dd = dir([verdn '/*.bam']);
        if length(dd)==1
          bamn = [verdn '/' dd(1).name];
        end
      end
      if exist(bamn,'file')
        bamd = dir(bamn);
        P.bampath{idx,1} = bamn;
        P.bamsize(idx,1) = bamd.bytes;
        P.bamtime{idx,1} = bamd.date;
        P.bamtimenum(idx,1) = bamd.datenum;
      else
        P.bampath{idx,1} = '?';
        P.bamsize(idx,1) = -1;
        P.bamtime{idx,1} = '?';
        P.bamtimenum(idx,1) = -1;
      end
      finishedfn = [verdn '/finished.txt'];
      if exist(finishedfn,'file')
        P.finished{idx,1} = 'true';
      else
        P.finished{idx,1} = 'false';
      end
    end % next version
  end % next sample
 end % next project
end % next rootdir

for f=1:length(flds)
  P = setfield(P,flds{f},F(:,f));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POST-CRAWL

% infer build
P.build = repmat({'?'},slength(P),1);
idx18 = grep('Homo_sapiens_assembly19',P.reference_sequence,1); P.build(idx18) = repmat({'hg19'},length(idx18),1);
idx19 = grep('Homo_sapiens_assembly18',P.reference_sequence,1); P.build(idx19) = repmat({'hg18'},length(idx19),1);
idx = grep('Canis_lupus_familiaris',P.reference_sequence,1); P.build(idx) = repmat({'dog'},length(idx),1);
idx = grep('Escherichia_coli',P.reference_sequence,1); P.build(idx) = repmat({'<i>E.coli</i>'},length(idx),1);
idx = grep(',',P.reference_sequence,1); P.build(idx) = P.reference_sequence(idx);

% impose initiative for certain projects
I = load_struct([rikerdir '/projects.txt']);
idx = listmap(P.project,I.project);
idx2 = find(~isnan(idx));
P.initiative(idx2) = I.initiative(idx(idx2));
P.initiative = fillblanks(P.initiative,'?');

% manually collapse CLL Tumor+Normal initiatives
P.initiative = regexprep(P.initiative,'^CLL (Normal|Tumor) WGS - Illumina$','CLL WGS - Illumina');

% choose records to ignore
P.ignore = false(slength(P),1);
igI1 = load_struct([rikerdir '/non-cancer_initiatives_to_ignore.txt']);
igI2 = load_struct([rikerdir '/not-this-build_initiatives_to_ignore.txt']);
igI = concat_structs({igI1,igI2});
P.ignore(ismember(P.initiative,igI.initiative)) = true;
igS = load_struct([rikerdir '/samples_to_ignore.txt']);
P.ignore(ismember(P.sample,igS.sample)) = true;

% unify certain initiative groups
ii = {'CLL WGS','CLL.*Hybrid Selection','Colon Cancer WGS','Esophageal Cancer.*WGS','Head.*Neck.*WGS','TCGA Lung Squamous WGS',...
   'Breast Cancer.*WGS','Breast Cancer.*Whole Exome.*SIGMA'};
for i=1:length(ii)
  idx = grep(ii{i},P.initiative,1);
  idx = idx(~P.ignore(idx));
  if ~isempty(idx)
    tmp = {};
    for i=1:length(idx), tmp = [tmp; split(P.initiative{idx(i)},',')]; end
    tmp = concat(unique(tmp),',');
    P.initiative(idx) = repmat({tmp},length(idx),1);
  end
end

% infer tech
P.tech = repmat({'capture'},slength(P),1);
idx = grep('WGS',P.analysis_type,1); P.tech(idx) = repmat({'WGS'},length(idx),1);
idx = grep('WGS',P.initiative,1); P.tech(idx) = repmat({'WGS'},length(idx),1);

% for BAMs that were created before the implementation of "finished.txt" (i.e. before 2009-07-16)
% we can impute finished status = true
idx = find(P.bamtimenum<734004.66635416669305413961 & P.bamtimenum>-1);
P.finished(idx) = repmat({'true'},length(idx),1);

% parse date out of timestamp
P.bamdate = regexprep(P.bamtime,'^(.*) \d\d:\d\d:\d\d$','$1');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finally,
% for samples in projects that have informative names (e.g. TCGA), 'guess' the individual and T/N

P.guess_individual = repmat({'?'},slength(P),1); P.guess_tn = P.guess_individual;

% TCGA
tmp = parse(P.sample,'^(TCGA-..-....)-([01])(.*)',{'indiv','tn','suffix'});
idx = find(~cellfun('isempty',tmp.indiv) & ~cellfun('isempty',tmp.tn));
P.guess_individual(idx) = tmp.indiv(idx); P.guess_tn(idx) = map_across(tmp.tn(idx),{'0';'1'},{'tumor';'normal'});
map=[];  map.name={'ovarian';'gbm';'aml'}; map.remove={'^TCGA';'^TCGA';'^TCGA-AB'};map.prefix={'OV';'GBM';'AML'};
for i=1:slength(map)
  idx2 = idx(grepi(map.name{i},P.initiative(idx),1));
  P.guess_individual(idx2) = regexprep(P.guess_individual(idx2),map.remove{i},map.prefix{i});
end
% for TCGA capture samples containing a "D" notation (indicating that the sample derives from a Broad-performed WGA),
% append "-BW" to the imputed individual name.
idx = intersect(grepi('capture|hybrid',P.initiative,1),grepi('^..-..D',tmp.suffix,1));
P.guess_individual(idx) = regexprep(P.guess_individual(idx),'(.*)','$1-BW');

% samples ending in "T" or "N"
tmp = parse(P.sample,'(.+[^_])_*(T|t|N|n)$',{'indiv','tn'});
idx = find(~cellfun('isempty',tmp.indiv) & ~cellfun('isempty',tmp.tn));
P.guess_individual(idx) = tmp.indiv(idx); P.guess_tn(idx) = map_across(lower(tmp.tn(idx)),{'t';'n'},{'tumor';'normal'});

% samples ending in "Tumor" or "Normal" or "Blood"
tmp = parse(P.sample,'^(.+)[-_]*(Tumor|tumor|Normal|normal|Blood|blood)$',{'indiv','tn'});
tmp.tn = regexprep(tmp.tn,'(Blood|blood)','normal');
idx = find(~cellfun('isempty',tmp.indiv) & ~cellfun('isempty',tmp.tn));
P.guess_individual(idx) = tmp.indiv(idx); P.guess_tn(idx) = lower(tmp.tn(idx));

% samples beginning in "Tumor" or "Normal" or "Blood"
tmp = parse(P.sample,'^(Tumor|tumor|Normal|normal|Blood|blood)[-_]*(.+)$',{'tn','indiv'});
tmp.tn = regexprep(tmp.tn,'(Blood|blood)','normal');
idx = find(~cellfun('isempty',tmp.indiv) & ~cellfun('isempty',tmp.tn));
P.guess_individual(idx) = tmp.indiv(idx); P.guess_tn(idx) = lower(tmp.tn(idx));






















