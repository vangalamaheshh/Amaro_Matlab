function savestruct(S,filename)
%SAVESTRUCT save a struct array in a tab-delimited file
%
%   savestruct(S,FILENAME)
%
% Where S is the structure to save and FILENAME is the path to the
% output file. The column headers will be the names of the fields in S.

%! TODO test inputs, pass filter in
fields = fieldnames(S);
nf = length(fields);
cols = cell(1,nf);
for f=1:nf
    if isnumeric(S(1).(fields{f}))
        cols{f} = {fields{f},[S.(fields{f})]};
    else
        cols{f} = {fields{f},{S.(fields{f})}};
    end
end
write_filtered_tabcols(filename,[],cols{:});