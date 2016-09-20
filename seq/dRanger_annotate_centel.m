function X = dRanger_annotate_centel(X,margin,telsize)
% annotate for centromere/telomere overlap

cen = load_cen;
len = load_chrlen;
if ~exist('telsize','var'), telsize = 25000; end
if ~exist('margin','var'), margin = 1000000; end

nx = slength(X);
X.class1 = cell(nx,1);
X.class2 = cell(nx,1);
for i=1:nx
  if X.pos1(i)<=telsize+margin, X.class1{i} = 'tel';
  elseif X.pos1(i)>=len(X.chr1(i))-telsize-margin, X.class1{i} = 'tel';
  elseif X.pos1(i)>=cen(X.chr1(i),1)-margin &...
    X.pos1(i)<=cen(X.chr1(i),2)+margin, X.class1{i} = 'cen';
  else X.class1{i} = '-';
  end
  if X.pos2(i)<=telsize+margin, X.class2{i} = 'tel';
  elseif X.pos2(i)>=len(X.chr2(i))-telsize-margin, X.class2{i} = 'tel';
  elseif X.pos2(i)>=cen(X.chr2(i),1)-margin &...
    X.pos2(i)<=cen(X.chr2(i),2)+margin, X.class2{i} = 'cen';
  else X.class2{i} = '-';
  end
  if X.pos1(i)<=telsize, X.class1{i} = 'TEL';
  elseif X.pos1(i)>=len(X.chr1(i))-telsize, X.class1{i} = 'TEL';
  elseif X.pos1(i)>=cen(X.chr1(i),1) &...
    X.pos1(i)<=cen(X.chr1(i),2), X.class1{i} = 'CEN';
  end
  if X.pos2(i)<=telsize, X.class2{i} = 'TEL';
  elseif X.pos2(i)>=len(X.chr2(i))-telsize, X.class2{i} = 'TEL';
  elseif X.pos2(i)>=cen(X.chr2(i),1) &...
    X.pos2(i)<=cen(X.chr2(i),2), X.class2{i} = 'CEN';
  end
end
