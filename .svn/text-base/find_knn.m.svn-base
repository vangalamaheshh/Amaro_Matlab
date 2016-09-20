function [S,D]=find_knn(v,k,operation,frac,give_dist,is_dist,dist_type)
if nargin == 2
  operation = 'and'; 
  frac=0.1;
  give_dist=0;
  is_dist=0;
end

if ~strcmp(operation,'and')
  give_dist=0;
end

N=size(v,1);
S=sparse(N,N);
D=sparse(N,N);

%nv=sum(v.*v,2);
chunk=1000;
if ~is_dist
  for i=1:N
      if mod(i,100)==0
          disp(i);
      end
    j=mod(i-1,chunk)+1;
    if (j==1)
        if exist('dist_type','var')  
            dsts=dist(v(i:min(i+chunk-1,N),:),v,dist_type);
        else
            dsts=dist(v(i:min(i+chunk-1,N),:),v);
        end
    end
    dst=dsts(j,:);
    md=mean(dst);
    id2=find(dst<frac*md);
    if length(id2)<k+2
      %      warning(['not enough in id2 in ' num2str(i) ]);
      id2=1:N;
    end
    [sd,idx]=sort(dst(id2));
    S(i,id2(idx(1:(k+1))))=1;
    if (give_dist)
      ksd=sd(1:(k+1))';
      ksd(ksd<0)=0;
      D(i,id2(idx(1:(k+1))))=sqrt(ksd);
    end
    %   if mod(i,20) == 0
    %      i
    %   end
  end
else
  for i=1:N
    if exist('dist_type','var')  
        dst=dist(v(i,:),v,dist_type);
    else
        dst=dist(v(i,:),v);
    end 
    md=mean(dst);
    id2=find(dst<frac*md);
    if length(id2)<k+2
      %      warning(['not enough in id2 in ' num2str(i) ]);
      id2=1:N;
    end
    [sd,idx]=sort(dst(id2));
    S(i,id2(idx(1:(k+1))))=1;
    if (give_dist)
      ksd=sd(1:(k+1))';
      ksd(ksd<0)=0;
      D(i,id2(idx(1:(k+1))))=sqrt(ksd);
    end
    %   if mod(i,20) == 0
    %      i
    %   end
  end
  
end


switch lower(operation)
 case 'and'
  S=S.*S';
  D=D.*S;
 case 'or'
  S=spones(S+S');
 case 'none'
 otherwise
end



