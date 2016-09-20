function Dstruct = cn_to_Dstruct(cn_file)
% Dstruct = cn_to_Dstruct(cn_file)
% Opens a CN file and reads all data, putting it into a Dstruct object
%
% by Michael J.T. O'Kelly, 080415

%% Check inputs

varlist1 = {'cn_file'};

defaults = {'ERR'};

required = [1];

for idx = 1:length(varlist1)
    if required(idx) && (~exist(varlist1{idx},'var') || eval(['isempty(' varlist1{idx} ')']))
        error('Required input %s undefined.',varlist1{idx})
    elseif ~exist(varlist1{idx},'var') || eval(['isempty(' varlist1{idx} ')'])
        eval([varlist1{idx} '=' defaults{idx}])
    end
end

fid = fopen(cn_file);
header_line = fgetl(fid);
headers = textscan(header_line, '%s'); % Finds all tab-delimited headers
num_headers = length(headers{1});
format_string = ['%q%q%d' repmat('%f', 1, num_headers-3)];
columns = textscan(fid,format_string,'ReturnOnError',0, 'Headerlines', 0, 'TreatAsEmpty', 'NA');
fclose(fid);

columns{2} = regexprep(columns{2},'[Xx]','23');
columns{2} = regexprep(columns{2},'[Yy]','24');
columns{2} = cell2mat(cellfun(@str2num,columns{2},'UniformOutput',0));

Dstruct.chrn = columns{2};
Dstruct.pos = double(Dstruct.chrn).*1e11 + double(columns{3});  %expand position by chrn to allow sorting and matching
[Dstruct.pos,idx] = sort(Dstruct.pos);
Dstruct.chrn = Dstruct.chrn(idx);
Dstruct.marker = columns{1}(idx);
Dstruct.sdesc = headers{1}(4:num_headers)';
Dstruct.dat = cell2mat({columns{4:num_headers}});
if min(Dstruct.dat)>=0  % if not already log converted
	Dstruct.dat = log2(Dstruct.dat)-1;
end

clear columns;




