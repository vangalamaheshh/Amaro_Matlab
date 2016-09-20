function structable = read_R_table(fname,dlm,nantext)
% READ_R_TABLE read an R-style table as a Matlab struct
%
%  STRUCTABLE = read_R_table(FNAME,DLM,NANTEXT)
%
% Read an R-style delimited text file FNAME containg a data table with a
% header into a MATLAB struct, with the field names formed from the column
% headers. DLM specifies the delimiter, default tab=char(9). NANTEXT 
% specifies the token used for missing data, default 'NA'. Missing data
% will be signified by NaN in the struct array data.

% default tab-delimited table
if ~exist('dlm','var') || isempty(dlm)
    dlm = char(9);
end
% 
if ~exist('nantext','var')
    nantext = 'NA';
end

nantext = ['^',nantext,'$'];

cpp = read_dlm_file(fname,dlm);
comment_lines = cellfun(@(x) ~isempty(regexp(x{1},'^#','start','once')),cpp);
if any(comment_lines)
    cpp{comment_lines} = [];
end
cpp = vertcat(cpp{:});

if isempty(cpp)
    structable = [];
else

    % convert header into legitimate structure field name
    hdr = regexprep(cpp(1,:),'\s+','_');
    hdr = regexprep(hdr,'[^A-Za-z0-9_]+','');
    hdr = regexprep(hdr,'^([0-9]+)','X$1');

    ncols = length(hdr);
    nrows = size(cpp,1)-1;
    structcell = cell(2,ncols);
    structcell(1,:) = hdr;
    for col = 1:ncols
        trynum = str2num(char(regexprep(cpp(2:end,col),nantext,'NaN')));
        if length(trynum) == nrows
            structcell{2,col} = num2cell(trynum);
        else
            structcell{2,col} = cpp(2:end,col);
        end
    end
    structable = struct(structcell{:});
end
