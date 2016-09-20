function [newSI,Dstructcell,IDXsid] = write_dloc_to_sifile(Dstructfiles,SIfile,SIfilenew)
%WRITE_DLOC_TO_SIFILE Tell the sample info file and the SI struct where to
%find the data structures.
%
%   [NEWSI,DSTRUCTCELL] = WRITE_DLOC_TO_SIFILE(DSTRUCTFILE,SIFILE) returns
%   a cell array of data structures (DSTRUCTCELL) and updates the sample
%   info file and SI struct.  The updated sample info file and structure
%   NEWSI contains new fields that give the dat file of each array and the
%   index of the array's D struct in Dstructcell.  DSTRUCTFILE can be a
%   single file or a cell array of strings listing multiple data files.
%   SIfile must be a string referencing a single file.
%
%       History:
%
%           * 07 Sept 24 - Written by Jen Dobson (jdobson@broad.mit.edu)

%Load SIfile and make .dfile field
SI = read_sample_info_file(SIfile);
[SI.dfile] = deal('');
[SI.dnum] = deal([]);
IDXsid = [];

if ~iscell(Dstructfiles)
    Dstructfiles = {Dstructfiles};
    numDfiles = 1;
else
    numDfiles = length(Dstructfiles);
end

%Load and set D equal to the struct in Dstructfile
for k = 1:numDfiles
    dfile = Dstructfiles{k};
    load(dfile);
    dd = whos('-file',dfile);
    D = eval(dd.name);
    allDs{k} = D;
    [dum,match1,match2] = match_string_sets(D.sdesc,{SI.array});
    [SI(match2).dfile] = deal(dfile);
    [SI(match2).dnum] = deal(k);
    IDXsid = [IDXsid match2];
end

    


    
    
    
    
%% Write to new SI file


if ~exist('SIfilenew','var')
    temp = textscan(SIfile,'%[^.]');
    SIfilenew = [temp{1}{1} '_modified' datestr(now,'yymmdd') '.txt'];
end


siwrite(SIfilenew,SI);

newSI = SI;
Dstructcell = allDs;
