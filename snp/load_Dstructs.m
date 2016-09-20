function [AL,Dstructcell] = load_Dstructs(Dstructfiles,ALfile,AL,conserve_memory,...
    hdf5dir,cnvfile,byposition,log2_input)
%LOAD_DSTRUCTS - Load D structs from DSTRUCTFILES and convert to datastruct objects if necessary.
%
%   [NEWAL,DSTRUCTCELL] = load_Dstructs(DSTRUCTFILE,ALFILE,AL)
%   returns a cell array of data structures (DSTRUCTCELL) and updates the
%   AL struct.  The updated array list file and structure NEWAL contains 
%   new fields that give the dat file of each
%   array and the dat file index in ...  DSTRUCTFILE can be a single file
%   or a cell array of strings listing multiple data files. ALfile must be
%   a string referencing a single file.
%
%   The returned AL structure is a Nx1 struct array with N=number of
%   arrays(samples).  For each array, possible fields are:
%               array -- gives array name (REQUIRED)
%               dfile -- gives the name of the data file (relative to base
%               directory) where the .dat for that array is located.
%               dnum -- gives the index of the data file in
%               (if the data is located in Dstructfiles{n}, AL(k).dnum = n;

%! TODO the structure of this code is highly redundant and forces one to
%! make changes in two places (singular and plural D). There should be a
%! process_one_D subfunction.

%
%       History:
%
%           * 07 Sept 24 - Written by Jen Dobson (jdobson@broad.mit.edu)
%
%           * 5 Dec 07 -- Added option to convert to datastruct object.
%           (and renamed m-file)  (jdobson)
%
%           * 6 Dec 07 -- Removed write of new AL file.  (jdobson)
%           
%           * 4 March 11 -- added optional input byposition for call to
%           remove_cnv.m (barbarat)


%Load ALfile and make .dfile field
if ~exist('AL','var')
    AL = read_array_list_file(ALfile);
end
% defaults for optional variables
if ~exist('byposition','var')
    byposition = [];
end
if ~exist('log2_input','var') || isempty(log2_input)
    log2_input = 0;
end

totalarrays = 0;

if ~iscell(Dstructfiles)
    Dstructfiles = {Dstructfiles};
    numDfiles = 1;
else
    numDfiles = length(Dstructfiles);
end

%Load and set D equal to the struct in Dstructfile
allDs = {};
filecount = 0;
for k = 1:numDfiles  %loop over the Dfiles
    dfile = Dstructfiles{k};
    D = load_D(dfile);
    Dsize = size(D,2);
    if isempty(strmatch(filesep,dfile))
        cwd = pwd;
        dfile = [cwd filesep regexprep(dfile,[ regexp_filesep '.+' regexp_filesep],'')];  %make dfile into full path name for hdf5 files
    end
    
    if iscell(D)
        %loop over the dstruct cells in the loaded D
        for jj = 1:Dsize  %if the loaded data structure is a cell array of structs
             

            [dum,match1,match2] = match_string_sets(D{jj}.sdesc,{AL.array});
            
            if ~isempty(match1)
                
                [AL(match2).dfile] = deal(dfile);
                [AL(match2).dnum] = deal(k);
                D{jj} = reorder_D_cols(D{jj},match1);
                
                if exist('cnvfile','var') && ~isempty(cnvfile)
                    D{jj} = remove_cnv(D{jj},cnvfile,byposition);
                end

                D{jj} = order_by_pos(D{jj});
                
                if log2_input
                    % convert log2 input to linear space
                    D{jj}.dat = 2.^(1+D{jj}.dat);
                end
                
                filecount = filecount + 1;
                %CONSERVE MEMORY OPTIONS

                if conserve_memory
                    if ~strcmp(class(D{jj}.dat),'single')
                        D{jj}.dat = single(D{jj}.dat);
                    end
                    if ~strcmp(class(D{jj}.affy_calls),'uint8')
                        D{jj}.affy_calls(isnan(D{jj}.affy_calls)) = 4;
                        D{jj}.affy_calls = uint8(D{jj}.affy_calls);
                    end
                 
                    D{jj} = rmfield_if_exists(D{jj},'orig');
                    D{jj} = rmfield_if_exists(D{jj},'origidx');
                    D{jj} = rmfield_if_exists(D{jj},'gorigidx');
                end


                if conserve_memory==2 && ~isa(D{jj},'datastruct')
                    verbose('Converting structure from %s to datastruct object',30,dfile);
                    D{jj} = datastruct(D{jj});
                    dflds = intersect({'dat','affy_calls'},fieldnames(D{jj}));
                    dfile = ['D' num2str(filecount)];
                    D{jj} = convert_to_diskfield(D{jj},dflds,strcat(hdf5dir,dfile,'_',dflds,'.h5'),dflds);

                elseif conserve_memory==2 && isa(D{jj},'datastruct')
                    dflds = setdiff({'dat','affy_calls'},diskfieldnames(D{jj}));
                    dfile = ['D' num2str(filecount)];
                    D{jj} = convert_to_diskfield(Dstructs{k},dflds,strcat(hdf5dir,dfile,'_',dflds,'.h5'),dflds);

                end
                

                allDs{end+1} = D{jj};  %put arrays in same order as in AL, alD
                D{jj} = [];
                totalarrays = totalarrays + length(allDs{end}.sdesc);
                
            end
            
        end
    else
        % single file case
                [dum,match1,match2] = match_string_sets(D.sdesc,{AL.array});
                
                if ~isempty(match1)
                    [AL(match2).dfile] = deal(dfile);
                    [AL(match2).dnum] = deal(k);

                    D = reorder_D_cols(D,match1);  %put arrays in same order as in AL, alD
                    if exist('cnvfile','var') && ~isempty(cnvfile)
                        D = remove_cnv(D,cnvfile,byposition);
                    end
                    D = order_by_pos(D);
                    if log2_input
                        % convert log2 input to linear space
                        D.dat = 2.^(1+D.dat);
                    end
                    
                    filecount = filecount + 1;
                    %CONSERVE MEMORY OPTIONS

                    if conserve_memory
                        if ~strcmp(class(D.dat),'single')
                            D.dat = single(D.dat);
                        end
                        if isfield(D,'affy_calls') && ~strcmp(class(D.affy_calls),'uint8')
                            D.affy_calls(isnan(D.affy_calls)) = 4;
                            D.affy_calls = uint8(D.affy_calls);
                        end
                       
                        D = rmfield_if_exists(D,'orig');
                        D = rmfield_if_exists(D,'origidx');
                        D = rmfield_if_exists(D,'gorigidx');
                    end


                    if conserve_memory==2 && ~isa(D,'datastruct')
                        verbose('Converting structure from %s to datastruct object',30,dfile);
                        D = datastruct(D);
                        dflds = intersect({'dat','affy_calls'},fieldnames(D));
                        dfile = ['D' num2str(filecount)];
                        D = convert_to_diskfield(D,dflds,strcat(hdf5dir,dfile,'_',dflds,'.h5'),dflds);

                    elseif conserve_memory==2 && isa(D,'datastruct')
                        dflds = setdiff({'dat','affy_calls'},diskfieldnames(D));
                        dfile = ['D' num2str(filecount)];
                        D = convert_to_diskfield(Dstructs{k},dflds,strcat(hdf5dir,dfile,'_',dflds,'.h5'),dflds);

                    end



                    allDs{end+1} = D;




                    D = [];
                    totalarrays = totalarrays + length(allDs{end}.sdesc);
                end
    end
end

Dstructcell = allDs;

if totalarrays<length(AL)
    warning('Did not find a match to every sample in array list')%#ok
end


