function fhd = get_firehose_dir(samples)

if ~iscell(samples), samples = {samples}; flag=1; else flag = 0; end

fhbd = '/xchip/cga1/firehose_output/Individual';

for i=1:length(samples)
  name1 = upper(regexprep(samples{i},'/','-'));
  name1 = regexprep(name1,'(.*)-WU','WU-$1');
  name1 = regexprep(name1,'(.*)-BCM','BCM-$1');
  name1 = regexprep(name1,'-WGS$','/wgs');
  name1 = regexprep(name1,'-CAPTURE$','/capture');
  name1 = regexprep(name1,'LUNG-','LA-S0');
  name1 = regexprep(name1,'MEL-','ME');
  dir1 = [fhbd '/' name1];
  if ~exist(dir1,'dir')
    dir1a = regexprep(dir1,'MM-','MMRC');
    if ~exist(dir1a,'dir')
      error('Not found: %s or %s',dir1, dir1a);
    else
      dir1 = dir1a;
    end
  end
  fhd{i,1} = dir1;
end

fhd = regexprep(fhd,'/MM-(0421|0422|0425)/','/MM-$1-FIX/');

if flag, fhd = fhd{1}; end
