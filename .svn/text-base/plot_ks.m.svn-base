function [ks_p,val]=plot_ks(vec1,vec2,tail)

nan1=find(isnan(vec1));
nan2=find(isnan(vec2));

if ~isempty(nan1) | ~isempty(nan2)
    warning('Ignoring NaN values');
    vec1(isnan(vec1))=[];
    vec2(isnan(vec2))=[];
end

n1=length(vec1);
n2=length(vec2);
vec1=sort(vec1);
vec2=sort(vec2);
v=[vec1; vec2];
[u,i,j]=unique(v);

nu=length(u);
[u1,i1,j1]=unique([u; vec1]);
[u2,i2,j2]=unique([u; vec2]);

% figure(1);
% plot(u(j1((nu+1):end)),(1:n1)/n1,'r-'); hold on;
% plot(u(j2((nu+1):end)),(1:n2)/n2,'b-');

i1a=i1;
i1a(i1<=nu)=0;
i1a(i1>nu)=i1a(i1>nu)-nu;
c=i1a(1);
for i=2:length(i1a)
   i1a(i)=max(i1a(i),c);
   c=i1a(i);
end 
% 
i2a=i2;
i2a(i2<=nu)=0;
i2a(i2>nu)=i2a(i2>nu)-nu;
c=i2a(1);
for i=2:length(i2a)
   i2a(i)=max(i2a(i),c);
   c=i2a(i);
end 

i1b=i1a/n1;
i2b=i2a/n2;
plot(u,i1b,'r-','LineWidth',2); hold on;
plot(u,i2b,'b-','LineWidth',2);

if tail==0
    [val,xi]=max(abs(i1b-i2b));
elseif tail==1
    [val,xi]=max(i1b-i2b);
else    
    [val,xi]=max(i2b-i1b);
end

line([u(xi) u(xi)],[i1b(xi) i2b(xi)],'LineWidth',3,'Color','k');

[ks_h,ks_p,ks_stat]=kstest2(vec1,vec2,0.01,tail);
