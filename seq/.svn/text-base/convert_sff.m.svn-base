function X = convert_sff(A,P)
% X = convert_sff(A,P)
% 
% input is from sffread, array of structs with fields: Header (char), Sequence (uint8), Quality (uint8)
%
% output is struct of arrays with fields:
%    header (char): same as input
%    seq (char): same as input, but changed from uint8 to char
%    qual (char): same as input, but changed from uint8 to char by adding 33
%    trim (double): first base from the 3'-end that has quality<trim_qual_cutoff 
%
% Mike Lawrence 2010-04-20

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'trim_qual_cutoff',15);

na = length(A);
X=[];
X.header = cell(na,1); X.seq = cell(na,1); X.qual = cell(na,1); X.len = nan(na,1); X.trim = nan(na,1);
fprintf('Converting sff data... ');
for i=1:na, if ~mod(i,10000), fprintf('%d/%d ',i,na); end
  X.header{i} = A(i).Header;
  X.seq{i} = char(A(i).Sequence);
  q = A(i).Quality;
  X.qual{i} = char(q+33);
  X.len(i) = length(q);
  idx = find(q>=P.trim_qual_cutoff,1,'last');
  if isempty(idx), X.trim(i) = 0; else X.trim(i) = idx; end
end, fprintf('\n');

