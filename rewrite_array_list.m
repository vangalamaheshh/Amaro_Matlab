function status = rewrite_array_list(oldalf,newalf,translationfile)
%REWRITE_ARRAY_LIST rewrites an array list with different array names.
%
%STATUS = REWRITE_ARRAY_LIST(OLDALF,TRANSLATIONFILE)
%
%   OLDALF is a standard array list file.
%   TRANSLATIONFILE is a two column file with names from the array lsit in
%   the first column and the names they should be translated to in the
%   second column.
%
%       jdobson@broad.mit.edu
%
%---
%$Id$
%$Date: 2008-05-09 17:36:58 -0400 (Fri, 09 May 2008) $
%$LastChangedBy: jdobson $
%$Rev$

AL = read_array_list_file(oldalf);

fid = fopen(translationfile);
tf = textscan(fid,'%s%s');
fclose(fid);


[match,si,sj] = intersect({AL.array},tf{1});
if length(si)<length(AL)
    warning('Did not find translation match to all arrays')
end

fldnames = fieldnames(AL);
ALcell = struct2cell(AL);
ALcell(1,si) = tf{2}(sj);
AL = cell2struct(ALcell,fldnames,1);

infofilewrite(newalf,AL)