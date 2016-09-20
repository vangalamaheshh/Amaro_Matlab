function count = infofilewrite(filename,S)
%INFOFILEWRITE Write info file from structure.  
%   
%   COUNT = INFOFILEWRITE(FILENAME,S) writes an info structure S (usually an
%   array list file or sample info file), to file FILENAME placing fields
%   in columns with S field names as the column headers.  COUNT is the
%   number of bytes successfully written. Write operation overwrites
%   FILENAME if it already exists.


%% Open file to write
fid = fopen(filename,'w');


fidstruct = cell(length(S),1);
[fidstruct{:}] = deal(fid);

%% Write header row

fldnames = fieldnames(S);
line1format = [repmat('%s\t',1,length(fldnames)) '\n'];
fprintf(fid,line1format,fldnames{:});

%% Make string to pass to fprintf


%Make cell array listing the full fieldnames (i.e. 'S.ploidy')
structname = cell(length(fieldnames(S)),1);
[structname{:}] = deal('S.');   
flds = cellfun(@strcat,structname,fieldnames(S),'UniformOutput',0);

%Make cell array of double 2
catdim = repmat({2},length(flds),1);

%Add curly braces and commas
comma = repmat({'}'','},length(flds),1);
curbrace = repmat({'{'},length(flds),1);
commaflds = cellfun(@cat,catdim,curbrace,flds,comma,'UniformOutput',0);

%fpflds gives the arguments passed to fprintf in string form
fpflds = cat(2,commaflds{:});
fpflds = fpflds(1:end-1); %remove final comma



%% Make format string

form = cell(length(S),1);
[form{:}] = deal([repmat('%s\t',1,length(flds)) '\n']);


%% Write data

s = ['cellfun(@fprintf,fidstruct,form,' fpflds ')'];

count = eval(s);

fclose(fid);