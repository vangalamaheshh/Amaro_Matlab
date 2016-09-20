function D = change_data_pointer(D,cf,fields,newfiles,fname)
% CHANGE_DATA_POINTER saves datastruct disk fields to disk in HDF5 file and
% moves fielddata pointer to new dataset for reading and writing.
%
%   CHANGE_DATA_POINTER(D,CF,FIELDS,NEWFILE,fname)
%
%           CF -- boolean value to indicate whether diskfield should point
%           to new file.
%
%           FIELDS -- optional input to specify which fields to save to
%           disk
%
%
%           NEWFILES -- optional file names changes datastruct to point to
%           new file.  If no inputs are given, this function increments the
%           file number suffix.\
%
%           FNAME -- optional field to input the name of the file to look
%           for transfer from data.  (Use when calling from load_D2 to
%           prevent D.mat files that have been copied to new directories
%           from looking back to the original to get the transfer data.)
%
%
%

%       Revisions
%           29 Nov 07 -- Function created by Jen Dobson
%           (jdobson@broad.mit.edu)
%
%           24 Jan 07 -- Implemented dataset transfer with call to C executable.
%
%           31 Jan 07 -- added FNAME input.
%
%---
% $Id$
% $Date: 2008-08-13 13:35:27 -0400 (Wed, 13 Aug 2008) $
% $LastChangedBy: jdobson $
% $Rev$

if ~exist('fields','var')
    fields = diskfieldnames(D);
else
    if ~iscell(fields)
        if ischar(fields)
            fields = {fields};
        else
            error('FIELDS input must be a cell array of strings')
        end
    end
end


if isempty(fields)
    error('Could not find disk fields to save')
end

if ~exist('cf','var')
    cf = 0;
end



[intflds,dum,flidx] = intersect(fields,D.fieldnames);


fielddata = [D.fielddata{flidx}];


if exist('newfiles','var')
    if ischar(newfiles)
        newfiles = {newfiles};
    elseif ~iscell(newfiles)
        error('NEWFILES must be a cell array of strings')
    end

    if length(newfiles) ~= length(fields)
        error('Length of NEWFILES must be equal to length FIELDS')
    end

    newfiles = get_new_filenames(newfiles);
    
elseif ~exist('newfiles','var') && cf

    newfiles = get_new_filenames({fielddata.Filename});

else

    newfiles = {fielddata.Filename};
end


%% Write the data to new file location and change the fielddata in datastruct
%to reflect new location
for k = 1:length(flidx)

    idx = flidx(k);

    if exist(newfiles{k},'file')
        fileattrib(newfiles{k},'+w','u');
    end
        
    if exist('fname','var') && ~isempty(fname)
        fromfile = fname;
    else
        fromfile =  D.fielddata{idx}.Filename;
    end

    unixstr = ['~/CancerGenomeAnalysis/trunk/C/transh5tonew18 ' fromfile ' ' D.fielddata{idx}.Name ' ' newfiles{k} ];

    stat = unix(unixstr);    


    D.fielddata{idx} = get_dataset_info(newfiles{k},D.fielddata{idx}.Name);



end


