function Dtr=trim_data(D)

m=mean(D(:));
s=std(D(:));
Dtr=D;
Dtr(D<(m-3*s))=m-3*s;
Dtr(D>(m+3*s))=m+3*s;
