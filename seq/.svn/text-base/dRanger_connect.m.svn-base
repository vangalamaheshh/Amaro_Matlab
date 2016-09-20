function dRanger_connect(infile)

try

if ~exist(infile,'file'), error('Could not find file %s',infile); end

X = load_struct(infile);
X = reorder_struct(X,strcmp(X.normreads,'0'));
X = make_numeric(X,{'chr1','chr2','pos1','pos2','str1','str2','tumreads','normreads'});
nx = slength(X);
keyboard

margin = 2000;
for i=1:nx
  idx1 = find(X.chr1==X.chr1(i) & X.chr2==X.chr2(i));
  idx2 = idx1(find(abs(X.pos1(idx1)-X.pos1(i))<=margin & abs(X.pos2(idx1)-X.pos2(i))<=margin));
  if length(idx2)>1
    fprintf('R %d/%d\n',i,nx);
    for j=1:length(idx2), k=idx2(j);
      fprintf('  chr%d[%d] %d  -->  chr%d[%d] %d   (%d tumreads, %d normreads)\n',...
        X.chr1(k), X.str1(k), X.pos1(k),...
        X.chr2(k), X.str2(k), X.pos2(k),...
        X.tumreads(k), X.normreads(k));
    end
  end
end






catch me, excuse(me); end
