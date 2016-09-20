function cons = get_cons46(chr,pos,P)

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

fwb = ['/cga/tcga-gsc/home/lawrence/db/' P.build '/conservation46/all.fwb'];
demand_file(fwb);

cons = get_from_fwb(fwb,chr,pos,P.build);

cons(cons==200)=nan;








