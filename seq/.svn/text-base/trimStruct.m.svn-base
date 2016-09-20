function [X1]=trimStruct(X,k) 
f=fieldnames(X);
kmax=max(k);
X1=X;
if (isfield(X,'N'))
    if (X.N>0)
        N=X.N;
        if (length(N)>1)
            error(' field N should be element count ');
        end
    end
end
if ~exist('N','var')
    N=0;
    for n=1:length(f)
        f1=char(f(n));
        x=X.(f1);
        m=size(x);
        N=max(N,m(1));
    end
end

for n=1:length(f)
    f1=char(f(n));
    x=X.(f1);
    m=size(x); 
   % if ( (m(1)>0) && ~(strcmp(f1,'vlab')))
    if ( (m(1)==N) && ~(strcmp(f1,'vlab')) && ~(strcmp(f1,'N'))  && ~(strcmp(f1,'head')) && ~(strcmp(f1,'hd')) )
   % if ( (m(1)==N) && ~(strcmp(f1,'vlab'))  && ~(strcmp(f1,'head')) && ~(strcmp(f1,'hd')) )
       if (m(1)<kmax) 
            display(['ooops ' f1])
            continue;
            %return;
      end
      x1=x(k,:);
      X1.(f1)=x1;
    else
      X1.(f1)=x;    
    end
end
X1.N=length(k);
 
function test()
 
 area='~/Projects/MobileElement/test/sim2/Spanner/build'
 f=[area '/sim2_fused.20.retro.span']
 what='sim2'

 % load all retro fragments 
A=loadRetroSpan(f);
k=find(bitget(A.element,1)>0);
A1=trimSPanStruct(A,k);