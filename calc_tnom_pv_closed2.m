function [pv,v]=calc_tnom_pv_closed(n,p,s)
% v(p-s,n-s)/(n+p p)
A=p-s;
B=n-s;
C=p-n;
M=n+p;
if (A+B)==0
    v=1;
    pv=exp(-ln_nchoosek(n+p,p));
    return
end
b=ceil(M/(A+B));

vec1=[ zeros(1,b); ones(1,b)]; vec1=(vec1(:))';
vec2=[ ones(1,b); zeros(1,b)]; vec2=(vec2(:))';


v=0;
for w=1:b
    wrd=vec1(1:w);
    tv=((C-t(wrd,A,B)+M)/2);
    if (tv <= M) & (tv>=0)
%        v=v+(-1)^(w+1)*nchoosek(M,tv);
%        v=v+(-1)^(w+1)*nchoosek(M,tv)/nchoosek(n+p,p);
        v=v+(-1)^(w+1)*exp(ln_nchoosek(M,tv)-ln_nchoosek(n+p,p));
    end
    
    wrd=vec2(1:w);
    tv=((C-t(wrd,A,B)+M)/2);
    if (tv <= M) & (tv>=0)
%        v=v+(-1)^(w+1)*nchoosek(M,tv);
%        v=v+(-1)^(w+1)*nchoosek(M,tv)/nchoosek(n+p,p);   
        v=v+(-1)^(w+1)*exp(ln_nchoosek(M,tv)-ln_nchoosek(n+p,p));   
    end
end
%pv=v/nchoosek(n+p,p);
pv=v;
v=round(exp(log(pv)+ln_nchoosek(n+p,p)));


function t=t(wrd,A,B)
wrdab=zeros(1,length(wrd));
wrdab(find(wrd))=B;
wrdab(find(1-wrd))=A;
if wrd(end)==0
    t=2;
else
    t=-2;
end
t=t*sum(wrdab);

