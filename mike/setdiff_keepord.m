function x = setdiff_keepord(x,y);

x = unique_keepord(x);
y = unique(y);
x(ismember(x,y))=[];
