function [newAL,Dstructcell] = write_dloc_to_alfile(Dstructfiles,ALfile,AL)
%WRITE_DLOC_TO_ALFILE Tell the array list file and the AL struct where to
%find the data structures.
%
%   [NEWAL,DSTRUCTCELL,ALD] = WRITE_DLOC_TO_ALFILE(DSTRUCTFILE,ALFILE,AL)
%   returns a cell array of data structures (DSTRUCTCELL) and updates the
%   array list file and AL struct.  The updated array list file and
%   structure NEWAL contains new fields that give the dat file of each
%   array and the dat file index in ...  DSTRUCTFILE can be a single file
%   or a cell array of strings listing multiple data files. ALfile must be
%   a string referencing a single file.
%
%   The AL structure is a Nx1 struct array with N=number of
%   arrays(samples).  For each array, possible fields are:
%               array -- gives array name
%               SIidx -- gives  the index of that array when referencing
%               from SIidx  (SI(AL(k).SIidx).array = AL(k).array)
%               dfile -- gives the name of the data file (relative to base
%               directory) where the .dat for that array is located.
%               dnum -- gives the index of the data file in
%               (if the data is located in Dstructfiles{n}, AL(k).dnum = n;
%
%
%
%
%       History:
%
%           * 07 Sept 24 - Written by Jen Dobson (jdobson@broad.mit.edu)

%Load ALfile and make .dfile field
if ~exist('AL','var')
    AL = read_array_list_file(ALfile);
end

[AL.dfile] = deal('');
[AL.dnum] = deal([]);
alD = [];

totalarrays = 0;

if ~iscell(Dstructfiles)
    Dstructfiles = {Dstructfiles};
    numDfiles = 1;
else
    numDfiles = length(Dstructfiles);
end

%Load and set D equal to the struct in Dstructfile
allDs = {};
for k = 1:numDfiles  %loop over the Dfiles
    dfile = Dstructfiles{k};
    load(dfile);
    dd = whos('-file',dfile);
    D = eval(dd.name);
    Dsize = dd.size(2);
    clear dd
    if iscell(D)
        %loop over the dstruct cells in the loaded D
        for jj = 1:Dsize  %if the loaded data structure is a cell array of structs
            if ~strcmp(class(D{jj}.dat),'single')
                D{jj}.dat = single(D{jj}.dat);
            end
            if ~strcmp(class(D{jj}.affy_calls),'uint8')
                D{jj}.affy_calls(isnan(D{jj}.affy_calls)) = 4;
                D{jj}.affy_calls = uint8(D{jj}.affy_calls);
            end

            [dum,match1,match2] = match_string_sets(D{jj}.sdesc,{AL.array});
            [AL(match2).dfile] = deal(dfile);
            [AL(match2).dnum] = deal(k);
            allDs{end+1} = reorder_D_cols(D{jj},match1);  %put arrays in same order as in AL, alD
            D{jj} = [];
            totalarrays = totalarrays + length(allDs{end}.sdesc);
        end
    else
        if ~strcmp(class(D.dat),'single')
            D.dat = single(D.dat);
        end
        if ~strcmp(class(D.affy_calls),'uint8')
            D.affy_calls(isnan(D.affy_calls)) = 4;
            D.affy_calls = uint8(D.affy_calls);
        end

        [dum,match1,match2] = match_string_sets(D.sdesc,{AL.array});
        [AL(match2).dfile] = deal(dfile);
        [AL(match2).dnum] = deal(k);
        allDs{end+1} = reorder_D_cols(D,match1);  %put arrays in same order as in AL, alD
        D = [];
        totalarrays = totalarrays + length(allDs{end}.sdesc);

    end
end

       
    
%% Write to new AL file

warning('off','MATLAB:fprintf:InputForPercentSIsNotOfClassChar');
if exist('ALfile','var')
    temp = textscan(ALfile,'%[^.]');
    ALfilenew = [temp{1}{1} '_modified' datestr(now,'yymmdd') '.txt'];
    infofilewrite(ALfilenew,AL);
end

newAL = AL;
Dstructcell = allDs;

if totalarrays<length(AL)
    warning('Did not find a match to every sample in array list')%#ok
end


