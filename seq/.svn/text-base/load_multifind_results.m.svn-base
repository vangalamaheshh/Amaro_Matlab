function x = load_multifind_results(dirname)
% load_multifind_results(dirname)
%
% given multifind output directory <dirname>,
% retrieves the alignments.
%
% returns struct x with fields:
%   class (A,A2,B,C,D,E,F)
%   fle (000-999)
%   targ (which multifind target was matched)
%   id
%   chr1,start1,end1,strand1,editstring1
%   chr2,start2,end2,strand2,editstring2
%   nerr1,nerr2 -> number of errors (gaps+mismatches)
%   is_tumor, is_normal
%
% Mike Lawrence 2009-03-03

% RETRIEVE MULTIFIND RESULTS

d = dir([dirname '/*.out']);
x = cell(length(d),1);
fprintf('Loading multifind results: ');
for i=1:length(d)
  if ~mod(i,100), fprintf('%d/%d ',i,length(d)); end
  x{i} = load_struct([dirname '/' d(i).name],[repmat('%f',1,10) '%s%s'],0);
  x{i}.fle = repmat(str2double(d(i).name(1:3)),slength(x{i}),1);
  x{i}.class = repmat({d(i).name(4)},slength(x{i}),1);
end
fprintf('\n');
x = concat_structs(x);
x = rename_field(x,{'col1','col2','col3','col4','col5','col6','col7','col8','col9','col10','col11','col12'},...
  {'targ','id','chr1','strand1','start1','end1','chr2','strand2','start2','end2','editstring1','editstring2'});
x.nerr1 = editstring_to_n_errors(x.editstring1);
x.nerr2 = editstring_to_n_errors(x.editstring2);

% SPLIT CLASS E into:
%    CLASS E = from paired-end lane; other end could not be aligned.
%    CLASS F = from single-end lane

lanes = load_lanelist;
idx = find(~lanes.has_pair(x.fle));
x.class(idx) = repmat({'F'},length(idx),1);

% SPLIT CLASS A into:
%    CLASS A  = reads <10K apart, opposite strands (normal)
%    CLASS A2 = reads <10K apart, same strand (weird)

idx = find(strcmp(x.class,'A') & x.strand1==x.strand2);
x.class(idx) = repmat({'A2'},length(idx),1);

% MARK TUMOR vs. NORMAL

x.is_tumor = strcmp(lanes.TN(x.fle),'T');
x.is_normal = strcmp(lanes.TN(x.fle),'N');
