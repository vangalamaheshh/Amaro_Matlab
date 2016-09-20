function res=same(s1,s2)
%SAME Return true if two cell array of strings or string arrays are 
%identical, return false otherwsise.
%
%   RES = SAME(S1,S2) returns 1 if S1 and S2 are identical string arrays or
%   cell array of strings, returns 0 otherwise.
%
%   Example: s1 = {'dog','cat','fish','bird'};
%            s2 = {'cat','dog','fish','bird'};
%            s3 = {'dog','cat','fish','bird'};
%
%            res = same(s1,s2) returns res = 0
%            res = same(s1,s3) returns res = 1
%
if ~ischar(s1)
  s1=strvcat(s1);
end
if ~ischar(s2)
  s2=strvcat(s2);
end

res=(any(range(s1-s2))==0);
