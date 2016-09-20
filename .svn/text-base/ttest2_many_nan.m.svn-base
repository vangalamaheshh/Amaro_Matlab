function [p,ratio] = ttest2_many_nan(dat,cl1,cl2,tail,var0,no_p)
% if degrees of freedom is 0 then p=ratio=NaN (both classes have
% only <=1 (non NaN) sample
% if at least one class has 2 samples and the other 1, if there is
% a difference then p=1, ratio=+-Inf. If there is no difference
% p=ratio=NaN
% no_p = don't calc p

is_nan=isnan(dat);
if any(is_nan(:))
  has_nan=1;

  notnan=~is_nan;
  nx=sum(notnan(:,cl1),2);
  ny=sum(notnan(:,cl2),2);  
else
  has_nan=0;
  nx=repmat(length(cl1),size(dat,1),1);
  ny=repmat(length(cl2),size(dat,1),1);
end


X=dat(:,cl1);
Y=dat(:,cl2);
clear dat;

if has_nan
  difference = nanmean(X,2) - nanmean(Y,2);
else
  difference = mean(X,2) - mean(Y,2);
end

% unknown but equal variances
if (0) 
  dfx = nx - 1; % length(cl1) -1 
  dfy = ny - 1; % length(cl2) - 1;
  dfe  = dfx + dfy;
  if exist('var0','var') && ~isempty(var0)
    %  msx = dfx .* (nanvar(X')'+repmat(var0,size(X,1),1));
    %  msy = dfy .* (nanvar(Y')'+repmat(var0,size(Y,1),1));
    if has_nan
      vx=nanvar(X,0,2);
    else
      vx=var(X,0,2);
    end    
    vx(vx<var0)=var0;
    if has_nan
      vy=nanvar(Y,0,2);
    else
      vy=var(Y,0,2);
    end
    vy(vy<var0)=var0;  
    msx = dfx .* vx;
    msy = dfy .* vy;
  else
    if has_nan
      msx = dfx .* nanvar(X,0,2);
      msy = dfy .* nanvar(Y,0,2);
    else
      msx = dfx .* var(X,0,2);
      msy = dfy .* var(Y,0,2);
    end
  end
  
  pooleds    = sqrt((msx + msy + eps) .* (ones(size(dfx,1),1)./(dfx +1) + ones(size(dfy,1),1)./(dfy + 1)) ./ dfe);
  ratio = difference ./ pooleds;
else
  if exist('var0','var') && ~isempty(var0)
    %  msx = dfx .* (nanvar(X,0,2)+repmat(var0,size(X,1),1));
    %  msy = dfy .* (nanvar(Y,0,2)+repmat(var0,size(Y,1),1));
    if has_nan
      vx=nanvar(X,0,2);
    else
      vx=var(X,0,2);
    end
    vx(vx<var0)=var0;
    if has_nan
      vy=nanvar(Y,0,2);
    else
      vy=var(Y,0,2);
    end
    vy(vy<var0)=var0;  
  else
    if has_nan
      vx = nanvar(X,0,2);
      vy = nanvar(Y,0,2);
    else
      vx = var(X,0,2);
      vy = var(Y,0,2);
    end
  end
  s2xbar = vx ./ nx;
  s2ybar = vy ./ ny;
  if (nx<=1) | (ny<=1)
    dfe=NaN;
  else
    dfe = (s2xbar + s2ybar) .^2 ./ (s2xbar.^2 ./ (nx-1) + s2ybar.^2 ./ (ny-1));
  end
  se = sqrt(s2xbar + s2ybar);
  ratio = difference ./ se;
end

%pooleds    = sqrt((msx + msy) .* (ones(size(dfx,1),1)./(dfx + 1) + ones(size(dfy,1),1)./(dfy + 1)) ./ dfe);

%disp(['difference=',num2str(difference)]);
%disp(['pooleds=',num2str(pooleds)]);


% ratio(find(difference==0))=0;
% Find the p-value for the tail = 1 test.
if no_p
  p=[];
else
  p  = 1 - tcdf(ratio,dfe);
  
  if ~exist('tail','var')
    tail=0;
  end
  % Adjust the p-value for other null hypotheses.
  if (tail == 0)
    p = 2 * min(p, 1-p);
    %    if (nargout>2)
    %        spread = tinv(1 - alpha / 2,dfe) * pooleds;
    %        ci = [(difference - spread) (difference + spread)];
    %    end
  else
    if (tail == 1)
      %       if (nargout>2)
      %           spread = tinv(1 - alpha,dfe) * pooleds;
      %           ci = [(difference - spread), Inf];
      %       end
    else
      p = 1 - p;
      %       if (nargout>2)
      %           spread = tinv(1 - alpha,dfe) * pooleds;
      %           ci = [-Inf, (difference + spread)];
      %       end
    end
  end
end




