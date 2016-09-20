function [gene strand type] = mutx(R,chr,pos)
% mutx(R,m)
% R = load_refseq(build);
% Mike Lawrence 2009-08-07

build ='hg18'

nr = length(chr);
if length(pos)~=nr, error('chr and pos must be same length'); end

gene = cell(nr,1);
strand = cell(nr,1);
type = cell(nr,1);
for i=1:nr, if ~mod(i,10), fprintf('%d/%d ',i,nr); end
  A = find_mut_in_refseq(R,build,['chr' num2str(chr(i))],pos(i),pos(i),'a');
  if slength(A)==0 | isempty(A.gene{1})
    type{i} = 'IGR';
    gene{i} = '';
    strand{i} = '';
  else
    gene{i} = A.gene{1};
    strand{i} = A.strand{1};
    if isempty(A.proteinchange{1})
      type{i} = 'intron';
    else
      type{i} = 'coding';
    end
  end
end
