function [ma,mb,dab,fa,fb]=matchI(a0,b0,W,F,C)
%  [ma,mb,dxab,fa,fb]=match2(a,b,w,f,c)
% match two lists of (x1,x2) coordinates within a window W and reciprocal  overlap F
% C=1 requires b entirely inside a (or b=a). 
if(nargin<5)
    C=0;
end
haswidth=true;
if ~isfield(a0, 'x1')
    a.x1=a0; 
    haswidth=false;
else
    a=a0;
end
if ~isfield(b0, 'x1')
    b.x1=b0; 
    haswidth=false;
else
    b=b0;
end
% length of a
Na=length(a.x1);
if ~isfield(a, 'x2')
    a.x2=a.x1+1; 
    haswidth=false;
end
% length of a
Nb=length(b.x1);
if ~isfield(b, 'x2')
    b.x2=b.x1+1; 
    haswidth=false;
end
% initialize
dab=[];
ma=0*a.x1;
mb=0*b.x1;
fa=ma; 
fb=mb;


% add W to b  (matching slosh)
b.x2=round(b.x2+W);
b.x1=round(b.x1-W);

% bx1=int64(round(b.x1));
% bx2=int64(round(b.x1));

% calc lengths of events
a.L=a.x2-a.x1+1;
b.L=b.x2-b.x1+1;

Nmatch=0;
% loop over a:  map a to b
for n=1:Na
    % some overlap with b
    ax1=round(a.x1(n));
    ax2=round(a.x2(n));    
    o=(b.x1<=ax2)&(ax1<=b.x2);
    if ~any(o), continue; end
    ko=find(o);
    bo.x1=b.x1(ko);
    bo.x2=b.x2(ko);
    bo.L=b.L(ko);
    
    % [ax1 bo.x1 ax2 bo.x2]-1e10
    
    % a before b
    o1 = (ax1<=bo.x1)  &  (ax2<bo.x2);  
    %o1 = (ax1<bo.x1)  &  (ax2<bo.x2);  
    %o1 = (ax1<bo.x1)  &  (ax2<=bo.x2);  
    %o1 = (a.x1(n)<bo.x1)  &  (a.x2(n)<bo.x2);  
    f1  = 0*o1;
    if any(o1)
        k1=find(o1);
        f=double(a.x2(n)-bo.x1(k1)+1)./double(max(a.L(n),bo.L(k1)));
        f1(k1)=f;
    end
    
    % b inside a or b=a
    o2 = (ax1<=bo.x1)  &  (ax2>=bo.x2);  
    %f2=(b.x2(ko)-b.x1(ko))./max(a.L(n),b.L(ko));
    f2  = 0*o2;
    if any(o2)
        k2=find(o2);
        f=double(bo.x2(k2)-bo.x1(k2)+1)./double(max(a.L(n),bo.L(k2)));
        f2(k2)=f;
    end
    
    % a inside b
    o3 = (ax1>bo.x1)  &  (ax2<=bo.x2);  
    %f3=2*(a.x2(n)-a.x1(n))./(a.L(n)+b.L(ko));
    %f3=(a.x2(n)-a.x1(n))./max(a.L(n),b.L(ko));
    f3  = 0*o3;
    if any(o3)
        k3=find(o3);
        f=double(a.x2(n)-a.x1(n)+1)./double(max(a.L(n),bo.L(k3)));
        f3(k3)=f;
    end
   
    
    % b before a
    o4 = (ax2>bo.x2) & (ax1>bo.x1);    
    %o4 = (ax2>=bo.x2) & (ax1>bo.x1);    
    %f4=2*(b.x2(ko)-a.x1(n))./(a.L(n)+b.L(ko));
    %f4=(b.x2(ko)-a.x1(n))./max(a.L(n),b.L(ko));
    f4  = 0*o4;
    if any(o4)
        k4=find(o4);
        f=double(bo.x2(k4)-a.x1(n)+1)./double(max(a.L(n),bo.L(k4)));
        f4(k4)=f;
    end
  
    if (C==1)  % b inside a only (Contained)
        o1=0*o1;
        o3=0*o3;
        o4=0*o4;
        f2=1+0*f2;
    end
    
    % overlap fraction
    fo=(o1.*f1+o2.*f2+o3.*f3+o4.*f4);
    dx=(a.x1(n)+a.x2(n))/2-(b.x1(ko)+b.x2(ko))/2;
    if haswidth  % use reciprocal overlap for best match
        [f k1]=max(fo);
        k1=k1(f>F);
    else              % use minumim residual for best match
        [d k1]=min(dx);
        f =fo(k1);
    end
    if (~isempty(k1)) 
        k=ko(k1);
        Nmatch=Nmatch+1;
        ma(n)=k;
        mb(k)=n;
        dab(Nmatch)=(a.x1(n)+a.x2(n))/2-(b.x1(k)+b.x2(k))/2;
        fa(n)=f;
        fb(k)=f;
    end
end

function test()

fid = fopen('~/Projects/MobileElement/test/1000_alu_insertions.pos');
C = textscan(fid, '%d %*[^\n]','headerlines',1);  fclose(fid);
spike.x1=double(C{1});
 
  % positions of known alus
 fid = fopen('~/Projects/MobileElement/alus-genomewide/chr20_alus0.alu');
 C = textscan(fid, '%d %d %d %d%*[^\n]','headerlines',1);  fclose(fid);
 alu.x1=double(C{2});
 alu.x2=double(C{3});

 W=350; F=0
  [ma,mb,dxab,fa,fb]=matchI(spike,alu,W,F)
  
  xlim=0e7+[1 1e6];
  ks=find((spike.x1>xlim(1))&(spike.x1<xlim(2)));
  plot(spike.x1(ks)/1e6,spike.x1(ks)*0,'bs')
  ka=find((alu.x1>xlim(1))&(alu.x1<xlim(2)));
  x=[alu.x1(ka) alu.x2(ka)]/1e6;
  y=alu.x1(ka)*0; y=[y y];
  line(x',y','color','k','linestyle','-','linewidth',2)

  km=find((ma>0)&(spike.x1>xlim(1))&(spike.x1<xlim(2)));
  hold on;
  plot(spike.x1(km)/1e6,spike.x1(km)*0+0.1,'rs')
  hold off;  
  
  set(gca,'ylim',[-1 1]);
  
    function plotmatch2()
        plot(minmax([a.x1;a.x2;b.x1; b.x2]),[0 1+max([length(a.x1) length(b.x1)])],'w+',a.x1,1:length(a.x1),'bx',b.x1,0.5+(1:length(b.x1)),'rx')
        line([a.x1'; a.x2'], [1:length(a.x1);1:length(a.x1)],'color','b','linestyle','-')
        line([b.x1'; b.x2'], 0.5+[1:length(b.x1);1:length(b.x1)],'color','r','linestyle','-')
                 h  = findobj(gca,'Type','line');
        set(h,'linewidth',4)
        
  