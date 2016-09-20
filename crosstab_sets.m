function ct=crosstab_sets(varargin)
%CROSSTAB_SETS Generate cross-tabulation tables
%
%       CT = CROSSTAB_SETS(VARARGIN)
%
% Elements of varargin should be vectors of all the same length. CT is
% mxnxpx... cell array (table) where m is number of unique elements in
% varagin{1}, n is number of unique elements in varargin{2}, p is number of
% unique elements in varagin{3}, etc.  Each element (i,j,k,...) of cell
% array CT gives the indices of input vectors that correspond to the
% (ith,jth,kth,...) unique value of
% varargin{1},varargin{2},varargin{3},...
%
%   Example: crosstab_sets([1 3 2 1 2 3 1 1 3 3 2 2 1 3 2 2 2 2 1 1 3 1 3 2],...
%       [1 3 3 1 1 2 2 1 1 3 2 3 1 2 2 1 1 1 2 3 2 2 2 3])

sz=[];
N=length(varargin{1});
for i=1:length(varargin)
  s{i}=unique(varargin{i});
  sz(i)=length(s{i});
end

ct=cell(sz);
for i=1:prod(sz)
  tmp=1:N;
  idx=ind2sub_vec(sz,i);
  for j=1:length(varargin)
    tmp=intersect(tmp,find(varargin{j}==s{j}(idx(j))));
%    find(varargin{j}==s{j}(idx(j)))
  end
  ct{i}=tmp;
end

