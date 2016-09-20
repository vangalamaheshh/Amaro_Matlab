function X = getbamlist(srcdir)

if ~exist('srcdir','var'), error('Must supply srcdir'); end

mask ='*.aligned.duplicates_marked.bam';
d = dir([srcdir '/' mask]);
% remove broken links
discard = [];
for i=1:length(d)
  if isempty(d(i).bytes), discard=[discard;i]; end
end
d = d(setdiff(1:length(d),discard));

n = cell(length(d),1);for i=1:length(d), n{i}=d(i).name; end
X = parse(n,'TCGA-\d\d-(\d*)-(\d)..-...\.',{'patient','isnorm'});
X.isnorm = str2double(X.isnorm); p = unique(X.patient);
X.bamname = regexprep(n,'(.*)',[srcdir '/$1']);

p = unique(X.patient);
keep = p;
for i=1:length(p)
  tidx = find(strcmp(X.patient,p{i}) & ~X.isnorm);
  nidx = find(strcmp(X.patient,p{i}) & X.isnorm);
  if isempty(tidx) | isempty(nidx)
    fprintf('Skipping %s, which does not have tumor and normal BAMs.\n',p{i});
    keep = setdiff(keep,p{i});
    continue;
  end
  tbai = dir([regexprep(X.bamname{tidx},'bam$','') '*bai']);
  nbai = dir([regexprep(X.bamname{nidx},'bam$','') '*bai']);

  if isempty(tbai) || isempty(nbai)
    fprintf('Skipping %s, which does not have tumor and normal BAIs.\n',p{i});
    keep = setdiff(keep,p{i});
    continue;
  end
  if length(tbai)>1 || length(nbai)>1
    fprintf('Skipping %s, which has more than one tumor or normal BAI.\n',p{i});
    keep = setdiff(keep,p{i});
    continue;
  end
  fprintf('Found %s.\n',p{i});
  X.bainame{tidx,1} = [srcdir '/' tbai.name];
  X.bainame{nidx,1} = [srcdir '/' nbai.name];
end

X = reorder_struct(X,find(ismember(X.patient,keep)));
