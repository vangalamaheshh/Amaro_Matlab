function c = load_cons(chr,build)
%
% load_cons(chr,build)
%
% build: default = hg18
%

if ~exist('build','var'), build='hg18'; end

if isnumeric(chr)
  if chr==23, chr = 'chrX';
  elseif chr==24, chr = 'chrY';
  elseif chr==25, chr = 'chrM';
  else chr = ['chr' num2str(chr)];
  end
end

file = ['/xchip/tcga/gbm/analysis/lawrence/cons/' build '/' chr '.dat'];

fprintf('Loading conservation data from %s\n',file);

fid=fopen(file,'r');
c=fread(fid,inf,'uint8=>uint8');
fclose(fid);

