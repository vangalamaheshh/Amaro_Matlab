function [] = saveLegoStructToFile(name, XY, fid, isHeader)
% [] = saveLegoStructToFile(name, XY, fid)
%
%  Append a file with raw counts used to generate the lego plots.
%
% name -- name of the case
% 
% XY -- the lego struct 
%
% fid -- open file for writing
%
% isHeader = true or false whether to write out the header.
%

Nc=length(XY.n(:));


% Construct header
fmt1='%s';
if isHeader
    fprintf(fid,'%s','name');
end
for i=1:Nc
    fmt1=[fmt1 '\t%d'];
    if isHeader
        fprintf(fid,'\t%s',XY.cat{i});
    end
    
end
if isHeader
    fprintf(fid,'\n');
end
fmt1=[fmt1 '\n'];

% Write data
fprintf(fid,fmt1,name,XY.n(:));
