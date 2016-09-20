function [C_new,m1]=match_using_array_list_files(C,alf)
% MATCH_USING_ARRAY_LIST_FILES Decimate and return data structure with only 
% those samples listed in the array list file.
% 
%    [C_new,M1] = MATCH_USING_ARRAY_LIST_FILES(C,ALF) returns cell array of 
%    data structures C_new and cell array of index vectors M1 s.t. C{i}(m1{i}) 
%    gives the ith array list.  Input C is a cell array of data structures 
%    with required field .sdesc giving the sample names, and alf is a cell 
%    of array list filenames.
%
%---
% $Id$
% $Date$
% $LastChangedBy$
% $Rev$

if length(C)~=length(alf)
  error('C and alf lengths don''t match');
end

for i=1:length(C)
  al=read_dlm_file(alf{i});

  %al=al(1:2:end);

  [M,m1{i},m2{i}]=match_string_sets(C{i}.sdesc,cat(1,al{:}));

end

common_m=intersect_many(m2{:});



for i=1:length(C)
  C_new{i}=reorder_D_cols(C{i},m1{i}(find(ismember(m2{i},common_m))));
end

