function M = crosstab3(col1,col2,m,n,append1, append2)
  
% %as rows
%   if size(col1,2) ~= size(col2,2)
%     error('columns must be same length');
%   end
%   

if ~exist('m','var')
    m = max(col1);
end

if ~exist('n','var')
    n = max(col2);
end

if ~exist('append1','var')
    append1 = reshape(repmat(1:m,n,1),1,m*n);
end

if ~exist('append2','var')
    append2 = repmat(1:n,1,m);
end


 col1 = [col1 append1];
 col2 = [col2 append2];

idx = m*(col1-1)+col2;
idx = sort(idx);
M = diff([0 find([diff(idx) 1])]);
M = reshape(M,m,n)' - ones(m,n);





% % 
% % col1 = [col1 reshape(repmat(1:m,n,1),1,m*n)];
% % col2 = [col2 repmat(1:n,1,m)];
% [col2,s2idx] = sort(col2);
% col1 = col1(s2idx);
% [col1,s1idx] = sort(col1);
% col2 = col2(s1idx);
% 
% bp1 = col1.*logical([1 diff(col1)]);
% bp2 = col2.*logical([1 diff(col2)]);
% M = diff([find(bp1 | bp2) length(col1)+1]);
% 
% M = reshape(M,m,n)' - ones(m,n);


%as cols

% 
%   if size(col1,2) ~= size(col2,2)
%     error('columns must be same length');
%   end
%   
%   m = max(col1);
%   n = max(col2);
% 
% col1 = [col1 ; reshape(repmat((1:m),n,1),m*n,1)];
% col2 = [col2 ; repmat((1:n)',m,1)];
% 
% [col2,s2idx] = sort(col2);
% [col1,s1idx] = sort(col1(s2idx));
%  col2 = col2(s1idx);
% bp1 = col1.*logical([1; diff(col1)]);
% bp2 = col2.*logical([1; diff(col2)]);
% M = diff([find(bp1' | bp2') length(col1)+1]);
% M = reshape(M,m,n)' - ones(m,m);
% 

