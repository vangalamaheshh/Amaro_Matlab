function X=read_cel_text_file(fname)

f=fopen(fname);
form='%s%f%f%f%f';
x=textscan(f,form,'bufSize',10000000,'delimiter','\t','emptyValue',NaN,'treatAsEmpty','NA');
X.marker=x{1};
X.x=x{2};
X.y=x{3};
X.val=x{4};
X.std=x{5};
fclose(f);
