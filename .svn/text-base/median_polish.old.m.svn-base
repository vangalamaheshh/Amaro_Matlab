function x=median_polish(x,start_with,max_iter)

if ~exist('start_with','var')
  start_with='rows';
end
if ~exist('max_iter','var')
  max_iter=10;
end

vec=[];
if is_col(start_with)
  iter=1;
  while iter<=max_iter && nnz(vec)>0
    x=x-repmat(median(x,1),size(x,1),1);
    vec=median(x,2);
    x=x-repmat(vec,1,size(x,2));
    iter=iter+1;
    disp(iter);
  end
else
  iter=1;
  while iter<=max_iter && nnz(vec)>0
    x=x-repmat(median(x,2),1,size(x,2));
    vec=median(x,2); %% FIX ME
    x=x-repmat(vec,size(x,2),1);
    iter=iter+1;
    disp(iter);
  end  
end
