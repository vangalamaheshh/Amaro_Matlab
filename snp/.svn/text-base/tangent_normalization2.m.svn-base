function [Y,orth_res,coef,delta_mat]=tangent_normalization(X,N,orth_res,ceil_val,ignore_idx,slow_loop,matched_data)
% tangent_normalization(X,N)
% X - samples you want to normalize
% N - Normals to use. If N is missing or empty then treat X as if were normals and 
%     run a leave-one-out normalization

coef=[];

if ~exist('orth_res','var')
  orth_res=[];
end

if ~exist('ignore_idx','var')
  ignore_idx=[];
end

if ~exist('ceil_val','var')
  ceil_val=[];
end

if ~exist('slow_loop','var')
  slow_loop=[];
end


if ~exist('N','var') || isempty(N)
  if exist('slow_loop','var') && ~isempty(slow_loop) && slow_loop
    Y=zeros(size(X));
    for i=1:size(X,2)
      disp(i);
      q=[];
      r=[];
      e=[];
      Y(:,i)=tangent_normalization(X(:,i),X(:,setdiff(1:size(X,2),i)),orth_res,ceil_val,ignore_idx,slow_loop);
    end
  else
    m=mean(X,2);
    y=X-repmat(m,1,size(X,2));
    
    if exist('ceil_val','var') && ~isempty(ceil_val)
      y1=ceil_val*tanh((1./ceil_val)*y);
    else
      y1=y;
    end
    
    if exist('ignore_idx','var') && ~isempty(ignore_idx)
      y1(ignore_idx,:)=[];
    end
    
    if ~exist('orth_res','var') || isempty(orth_res)   
      % q=orth(y1);  
      [q,r,e]=qr(y1,0);
      q=q(:,1:(end-1));
      r=r(1:(end-1),:);
    end
    
    b=y;
    
    if exist('ceil_val','var') && ~isempty(ceil_val) && (ceil_val)
      b1=ceil_val*tanh((1./ceil_val)*b);
    else
      b1=b;
    end
    
    if exist('ignore_idx','var') && ~isempty(ignore_idx)
      b1(ignore_idx,:)=[];
    end
    
    prj=q'*b1;
    % FIX ME for use with ignore_idx and ceil_val
    % add option to use matching X and N normalization
    for i=1:size(X,2)
      disp(i)
      p=prj(:,i);
      d=-p/(size(X,2)-1);
      P=prj-repmat(d,1,size(prj,2));
      Pnew=P(:,setdiff(1:size(X1,2),i));
      % qnew=orth(Pnew);
      [qnew,rnew,enew]=qr(Pnew,0);
      qnew=qnew(:,1:(end-1));
      rnew=rnew(1:(end-1),:);
      Ynew=P(:,i)-qnew*(qnew'*P(:,i));
      Y(:,i)=q*Ynew; %#ok
    end 
  end
elseif exist('matched_data','var') && ~isempty(matched_data) && (matched_data==1)
  Y=zeros(size(X));
  for i=1:size(X,2)
    disp(i);
    q=[];
    r=[];
    e=[];
    Y(:,i)=tangent_normalization(X(:,i),N(:,setdiff(1:size(X,2),i)),orth_res,ceil_val,ignore_idx,slow_loop);
  end    
  orth_res=[];
  coef=[];
  delta_mat=[];
elseif exist('matched_data','var') && ~isempty(matched_data) && (matched_data==2)
  Y=zeros(size(X));
  m_all=nanmean(N,2);
  y_all=[X N]-repmat(m_all,1,size(N,2));
  
  if ~exist('orth_res','var') || isempty(orth_res)
    if exist('ceil_val','var') && ~isempty(ceil_val)
      y1_all=ceil_val*tanh((1./ceil_val)*y_all);
    else
      y1_all=y_all;
    end
    
    if exist('ignore_idx','var') && ~isempty(ignore_idx)
      y1_all(ignore_idx,:)=[];
    end
    
    % q=orth(y1);  
    [q_all,r_all,e_all]=qr(y1_all,0);
    q_all=q_all(:,1:(end-1));
    r_all=r_all(1:(end-1),:);
  else
    q_all=orth_res.q;
    r_all=orth_res.r;
    e_all=orth_res.e;
  end
  
  x_prj=q_all'*X;
 %%  X_all; %%%%%%% HERE : FIX ME 
 
  for i=1:size(X,2)
    
    disp(i);
    q=[];
    r=[];
    e=[];
    Y(:,i)=tangent_normalization(X(:,i),N(:,setdiff(1:size(X,2),i)),orth_res,ceil_val,ignore_idx,slow_loop);
  end    

else
  m=nanmean(N,2);
  y=N-repmat(m,1,size(N,2));
  if ~exist('orth_res','var') || isempty(orth_res)

    if exist('ceil_val','var') && ~isempty(ceil_val)
      y1=ceil_val*tanh((1./ceil_val)*y);
    else
      y1=y;
    end
    
    if exist('ignore_idx','var') && ~isempty(ignore_idx)
      y1(ignore_idx,:)=[];
    end
    
    % q=orth(y1);  
    [q,r,e]=qr(y1,0);
    q=q(:,1:(end-1));
    r=r(1:(end-1),:);
  else
    q=orth_res.q;
    r=orth_res.r;
    e=orth_res.e;
  end
  
  b=X-repmat(m,1,size(X,2));
  
  if exist('ceil_val','var') && ~isempty(ceil_val) && (ceil_val)
    b1=ceil_val*tanh((1./ceil_val)*b);
  else
    b1=b;
  end

  if exist('ignore_idx','var') && ~isempty(ignore_idx)
    b1(ignore_idx,:)=[];
  end
  
  [tmp,inv_e]=sort(e);
  coef=pinv(r)*(q'*b1);
  coef=coef(inv_e,:);
  
%  Y=b-y(:,e)*pinv(r)*(q'*b1);
  Y=b-y*coef;
  delta_mat=y*coef+repmat(m,1,size(X,2));
  orth_res.q=q;
  orth_res.r=r;
  orth_res.e=e;
end
