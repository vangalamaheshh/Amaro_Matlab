function liftover_maf(vargin)
% three ways to run:
%    liftover_maf(inmaf,outmaf): hg18->hg19 assumed
%    liftover_maf(inmaf,outmaf,inbuild,outbuild)
%    liftover_maf(inmaf,outmaf,inbuild,outbuild,libdir): Firehose-deployed

if nargin<2 || nargin>5
  error('usage: liftover_maf(inmaf,outmaf,build1,build2,libdir)');end
end

inmaf  = varargin{1};
outmaf = varargin{2};

fprintf('Input MAF: %s\n', inmaf);
fprintf('Output MAF: %s\n', outmaf);

if nargin==2
  build1=18; 
  build2=19;
  fprintf('Assuming you want hg18->hg19\n');
elseif nargin==4
  build1 = varargin{3};
  build2 = varargin{4};
elseif nargin==5
  build1 = varargin{3};
  build2 = varargin{4};
  libdir = varargin{5};
end

fprintf('Input build: %s\n', build1);
fprintf('Output build: %s\n', build2);

if exist('libdir','var')
  fprintf('libdir: %s\n', libdir);
  % set java classpath to make accessible all jars in the libdir
  javaclasspath([libdir;direc([libdir '/*.jar'])]);
end

b1 = interpret_build(build1);
b2 = interpret_build(build2);
if b1==18 && b2==19
  new_NCBI_Build = '37';
elseif b1==19 && b2==18
  new_NCBI_Build = '36';
else
  error('Unsupported combination of builds\n');
end
hg1 = ['hg' num2str(b1)];
hg2 = ['hg' num2str(b2)];

% load starting maf
x = load_struct(inmaf);

% save old coordinates
flds = {'Chromosome','Start_position','End_position','Genome_Change'};
old_suffix = ['_' hg1];
for i=1:length(flds)
  if isfield(x,flds{i})
    x = setfield(x,[flds{i} old_suffix],getfield(x,flds{i}));
end,end

% get new coordinates
chr = x.Chromosome;
oldstart = x.Start_position;
oldend = x.End_position;
oldstartnum = str2double(oldstart);
oldendnum = str2double(oldend);
x.Start_position = liftover(chr,oldstartnum,build1,build2);
x.End_position = liftover(chr,oldendnum,build1,build2);
x.NCBI_Build = repmat({new_NCBI_Build},slength(x),1);

if slength(x)>0

  % find bad mappings
  bad = isnan(x.Start_position) | isnan(x.End_position);
  bad(oldstartnum>oldendnum) = true;
  bad(x.Start_position>x.End_position) = true;
  bad((oldendnum-oldstartnum)~=(x.End_position-x.Start_position)) = true;

  % also check for mismatched reference
  fprintf('Checking reference: ');
  for i=1:slength(x),if ~mod(i,100), fprintf('%d/%d ',i,slength(x)); end
    if bad(i), continue; end
    ref = x.Reference_Allele{i};
    if strcmp(ref,'-'), continue; end  % insertion
    oldref = genome_region(chr{i},oldstartnum(i),oldendnum(i),hg1);
    newref = genome_region(chr{i},x.Start_position(i),x.End_position(i),hg2);
    if ~strcmpi(oldref,ref), fprintf('ref %s oldref %s\n',ref,oldref); bad(i)=true; end
    if ~strcmpi(newref,ref), fprintf('ref %s newref %s\n',ref,newref); bad(i)=true; end
  end, fprintf('\n');

  % mark bad mapping as "Invalid"
  if any(bad)
    x.Validation_Status(bad) = repmat({'Invalid:failed_liftOver'},sum(bad),1);
    x.Chromosome(bad) = repmat({'Unknown'},sum(bad),1);
    x.Start_position(bad) = 0;
    x.End_position(bad) = 0;
  end

  % update Genome_Change (if field exists)
  if ~isfield(x,'Genome_Change')
    if isfield(x,'genomechange')
      x = rename_field(x,'genomechange','Genome_Change');
    else
      fprintf('(no GenomeChange field found)\n');
    end
  end
  if isfield(x,'Genome_Change')
    for i=1:slength(x)
      x.Genome_Change{i} = regexprep(x.Genome_Change{i},chr{i},x.Chromosome{i});
      x.Genome_Change{i} = regexprep(x.Genome_Change{i},oldstart{i},num2str(x.Start_position(i)));
      x.Genome_Change{i} = regexprep(x.Genome_Change{i},oldend{i},num2str(x.End_position(i)));
    end
  end
end

% save finished maf
save_struct(x,outmaf);

if exist('libdir','var')     % if running in Firehose-deployed mode, then:
  close all                  % close all figures so xvnc can terminate
end
