function s = get_splicedist(chr,pos,P)

if exist('P','var') && ischar(P)
  fprintf('Assuming %s is build\n',P);
  tmp=P;
  P=[];
  P.build = tmp;
  clear tmp
end

if ~exist('P','var'), P=[]; end

if ~isfield(P,'build')
  fprintf('Assuming hg19\n');
  P.build = 'hg19';
end

if ~isnumeric(chr), chr = convert_chr(chr); end
if ~isnumeric(pos), pos = str2double(pos); end

fwb = ['/cga/tcga-gsc/home/lawrence/db/' P.build '/splicedist/all.fwb'];
demand_file(fwb);

s = get_from_fwb(fwb,chr,pos,P.build);

s(s>2e9) = -(s(s>2e9)-2e9);


