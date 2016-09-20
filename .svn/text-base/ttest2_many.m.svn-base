function [p,ratio] = ttest2_many(dat,cl1,cl2,tail,var0)


dfx = length(cl1) - 1;
dfy = length(cl2) - 1;
dfe  = dfx + dfy;
X=dat(:,cl1);
Y=dat(:,cl2);
if nargin==5
  msx = dfx * (var(X')'+repmat(var0,size(X,1),1));
  msy = dfy * (var(Y')'+repmat(var0,size(Y,1),1));
else
  msx = dfx * var(X')';
  msy = dfy * var(Y')';
end

difference = mean(X,2) - mean(Y,2);
pooleds    = sqrt((msx + msy + eps) * (1/(dfx + 1) + 1/(dfy + 1)) / dfe);

ratio = difference ./ pooleds;

ratio(find(difference==0))=0;
% Find the p-value for the tail = 1 test.
p  = 1 - tcdf(ratio,dfe);

if nargin==3
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


