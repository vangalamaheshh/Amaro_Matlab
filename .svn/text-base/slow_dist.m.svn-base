function [d,res]=slow_dist(X,Y,dist_type,slow_data,transpose_X,transpose_Y)

res={};
if exist('transpose_X','var') && transpose_X
  xT=1;
  nX=size(X,2);
else
  xT=0;
  nX=size(X,1);
end

if exist('transpose_Y','var') && transpose_Y
  yT=1;
  nY=size(Y,2);
else
  yT=0;
  nY=size(Y,1);
end

d=zeros(nX,nY);

if slow_data==1
  if yT
    Y=Y';
  end
  for i=1:nX
    if xT
      d(i,:)=dist(X(:,i)',Y,dist_type,0,0);
    else
      d(i,:)=dist(X(i,:),Y,dist_type,0,0);
    end
    if mod(i,100)==0
      disp(i);
    end
  end
elseif slow_data==2
  if xT
    X=X';
  end
  for i=1:nY
    if yT
      d(:,i)=dist(X,Y(:,i)',dist_type,0,0);
    else
      d(:,i)=dist(X,Y(i,:),dist_type,0,0);
    end      
    if mod(i,100)==0
      disp(i);
    end
  end
else
  error('slow_data should be 1 or 2');
end

