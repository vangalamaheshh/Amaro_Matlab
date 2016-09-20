function [outputStruct,vals]=geosoftread(geotext)
% GEOSOFTREAD reads in Gene Expression Omnibus SOFT format data.
%
%   GEODATA = GEOSOFTREAD(FILE) reads in Gene Expression Omnibus (GEO) SOFT
%   format Sample (GSM), DataSet (GDS), and Platform (GPL) data from FILE
%   and creates a structure GEODATA, containing the following fields:
%       Scope
%       Accession
%       Header
%       ColumnDescriptions
%       ColumnNames
%       Data
%
%   GDS data will also have fields Identifier and IDRef.
%
%   FILE can also be a URL or a MATLAB character array that contains the
%   text of a GEO SOFT format file.
%
%   Example:
%
%       % Get a file from GEO and save it to a file.
%       geodata = getgeodata('GSM3258','TOFILE','GSM3258.txt')
%
%       % In subsequent MATLAB sessions you can use geosoftread to access the
%       % local copy from disk instead of accessing it from the GEO web site.
%       geodata = geosoftread('GSM3258.txt')
% 
%       % Read the expression GDS file for photosynthesis in proteobacteria.
%       gdsdata = geosoftread('GDS329.soft')
%
%   See also AFFYREAD, AGFEREAD, GALREAD, GETGEODATA, GPRREAD, SPTREAD.

% Copyright 2003-2008 The MathWorks, Inc.
% $Revision: 1.8.6.17 $   $Date: 2008/06/24 16:59:57 $

bioinfochecknargin(nargin,1,mfilename);

if ~ischar(geotext) && ~iscellstr(geotext)
    error('Bioinfo:geosoftread:InvalidStringInput',...
        'The input is not an array of characters or a cell array of strings.');
end

% If the input is a string of GenBank data then first character should be ^
isAstring = ischar(geotext) && ~isempty(geotext) && geotext(1) == '^';

bufSize = 2^16;
if ~isAstring
    if ~iscellstr(geotext)
        if (strfind(geotext, '://'))
            if (~usejava('jvm'))
                error('Bioinfo:geosoftread:NoJava','Reading from a URL requires Java.')
            end
            % must be a URL
            geotext = urlread(geotext);
            % clean up any &amp s
            geotext=strrep(geotext,'&amp;','&');
            geotext = textscan(geotext,'%s','delimiter','\n','Whitespace','','BufSize',bufSize);
            geotext = geotext{1};
        else
            if exist(geotext,'file')
                % it is a file
                fid = fopen(geotext,'rt');
                geotext = textscan(fid,'%s','delimiter','\n','Whitespace','','BufSize',bufSize);
                fclose(fid);
                geotext = geotext{1};
            elseif ~exist(geotext) %#ok
                error('Bioinfo:geosoftread:FileNotExist','The file %s does not exist.',geotext)
            end
        end
    else
        error('Bioinfo:geosoftread:InvalidSOFTFile','FILE is not a valid SOFT file.')
    end
else
    geotext = strread(geotext,'%s','delimiter','\n');
end


%line number
% ln = 1;

emptyLines = cellfun('isempty',geotext);
geotext(emptyLines) = [];
numLines = size(geotext,1);
minusTab = sprintf('-\t');
minusZeroTab = sprintf('-0\t');

% format as follows
% ^ start line with ID
% ! comment lines
% # column descriptions
% data values

colHeaderInfo = strncmp(geotext,'#',1);
firstColHeader = find(colHeaderInfo,1);
startLn = max(firstColHeader -1,1);

GEOType = strtrim(lower(strtok(geotext{1}(2:end),'=')));

% SAMPLE data
switch GEOType
    case {'sample','database','platform'}
        try
            if strcmp(GEOType,'database')
                dsLines = strmatch('^dataset',lower(geotext));
                isDB = true;
            else
                dsLines = 1;
                isDB = false;
            end
            [outputStruct.Scope, outputStruct.Accession] = strtok(geotext{dsLines(1)},'=');
            outputStruct.Scope = strtrim(outputStruct.Scope(2:end));
            outputStruct.Accession = strtrim(outputStruct.Accession(2:end));
            ln = startLn+1;
            % step over comment lines
            while geotext{ln}(1) == '!'
                geotext{ln}(1) = '';
                ln = ln+1;
            end
            outputStruct.Header.Type ='Gene Expression Omnibus';
            outputStruct.Header.Text = char(geotext(2:ln-1));

            colStart = ln;
            while geotext{ln}(1) == '#'
                geotext{ln}(1) = '';
                ln = ln+1;
            end

            outputStruct.ColumnDescriptions = geotext(colStart:ln-1);
            % get rid of any comment lines
            while geotext{ln}(1) == '!'
                geotext{ln}(1) = '';
                ln = ln+1;
            end
            % Attempt to guess the format of the columns -- numeric or text
            outputStruct.ColumnNames = strread(geotext{ln},'%s','delimiter','\t');
            numCols = numel(outputStruct.ColumnDescriptions);
            numRows = numLines-ln;
            geotext{ln+1} = strrep(geotext{ln+1},'Error','NaN');
            geotext{ln+1} = strrep(geotext{ln+1},'error','NaN');
            geotext{ln+1} = strrep(geotext{ln+1},'null','NaN');
            geotext{ln+1} = strrep(geotext{ln+1},minusTab,minusZeroTab);
            splitLine = strread(geotext{ln+1},'%s','delimiter','\t');
            numRead = numel(splitLine);
            isNumeric = [true(1,numRead) false(1,numCols-numRead)];
            for count = 1:numRead
                if isempty(splitLine{count}) || ~isempty(regexp(splitLine{count},'[^0-9.eE\-]','once')) && ~strcmp(splitLine{count},'NaN')
                    isNumeric(count) = false;
                end
            end
            if isDB
                isNumeric(2) = false;
            end
            badDataWarning = false;
            goodLines = true(numRows,1);
            if all(isNumeric)
                vals = zeros(numRows,numCols);
                for theLine = 1:numRows
                    % deal with comments or blank lines
                    if length(strtrim(geotext{theLine+ln})) < 1 || geotext{theLine+ln}(1) == '!' || geotext{theLine+ln}(1) == char(22)
                        goodLines(theLine) = false;
                        continue
                    end
                    try
                        vals(theLine,:) = strread(geotext{theLine+ln},'%f','delimiter','\t','emptyvalue',NaN)';
                    catch allExceptions %#ok<NASGU>
                        geotext{theLine+ln} = strrep(geotext{theLine+ln},'Error','NaN');
                        geotext{theLine+ln} = strrep(geotext{theLine+ln},'error','NaN');
                        geotext{theLine+ln} = strrep(geotext{theLine+ln},'null','NaN');
                        geotext{theLine+ln} = strrep(geotext{theLine+ln},'#DIV/0!','NaN');
                        geotext{theLine+ln} = strrep(geotext{theLine+ln},minusTab,minusZeroTab);
                        try
                            vals(theLine,:) = strread(geotext{theLine+ln},'%f','delimiter','\t','emptyvalue',NaN)';
                        catch moreExceptions %#ok<NASGU>
                            vals(theLine,:) = nan(1,numCols);
                            badDataWarning = true;
                        end
                    end
                end
                vals = vals(goodLines,:);
            else
                vals = cell(numRows,numCols);
                percents = repmat('%',1,numCols);
                Fs = repmat('s',1,numCols);
                Fs(isNumeric) = 'f';
                formatString = reshape([percents;Fs],1,2*numCols);
                nonNumericCols = find(~isNumeric);
                for theLine = 1:numRows

                    if isempty(geotext{theLine+ln})
                        continue
                    end
                    try

                        % deal with comments or blank lines
                        if length(strtrim(geotext{theLine+ln})) < 1 || geotext{theLine+ln}(1) == '!' || geotext{theLine+ln}(1) == char(22)
                            goodLines(theLine) = false;
                            continue
                        end
                        [vals{theLine,:}] = strread(geotext{theLine+ln},formatString,'delimiter','\t','emptyvalue',NaN);
                    catch allExceptions %#ok<NASGU>
                         geotext{theLine+ln} = strrep(geotext{theLine+ln},'Error','NaN');
                        geotext{theLine+ln} = strrep(geotext{theLine+ln},'error','NaN');
                        geotext{theLine+ln} = strrep(geotext{theLine+ln},'null','NaN');
                        geotext{theLine+ln} = strrep(geotext{theLine+ln},'NULL','NaN');
                        geotext{theLine+ln} = strrep(geotext{theLine+ln},'#DIV/0!','NaN');
                        geotext{theLine+ln} = strrep(geotext{theLine+ln},minusTab,minusZeroTab);
                        try
                            [vals{theLine,:}] = strread(geotext{theLine+ln},formatString,'delimiter','\t','emptyvalue',NaN);
                        catch moreExceptions %#ok<NASGU>
                            vals(theLine,isNumeric) = {NaN};
                            vals(theLine,~isNumeric) = {''};
                            badDataWarning = true;
                        end
                    end
                    for cols = nonNumericCols
                        vals{theLine,cols} = char(vals{theLine,cols});
                    end
                end
            end
            if any(~goodLines)
                vals = vals(goodLines,:);
            end
            if isDB
                if  strcmpi(outputStruct.ColumnNames{1},'ID_REF') &&...
                        strcmpi(outputStruct.ColumnNames{2},'IDENTIFIER')
                    outputStruct.IDRef = vals(:,1);
                    outputStruct.Identifier= vals(:,2);
                    outputStruct.ColumnNames(1:2) = [];
                    outputStruct.ColumnDescriptions(1:2) = [];
                    startCol = 3;
                end
            else
                startCol =1;
            end

            % see if we have a numeric array with spaces and try to
            % squeeze this into a matrix rather than a cell.
            for cols = startCol:numCols
                if ~isNumeric(cols)
                    tempCol = vals(:,cols);
                    strData = false;
                    for count = 1:numel(tempCol)
                        tempCol{count} = str2double(tempCol{count});
                        if isnan(tempCol{count}) && ~isempty(vals{count,cols})
                            strData = true;
                            break
                        end
                    end
                    if strData == false
                        vals(:,cols) = tempCol;
                        isNumeric(cols) = true;
                    end
                end
            end
            if all(isNumeric(startCol:end))
                vals = vals(:,startCol:end);
                if iscell(vals)
                    vals = cell2mat(vals);
                end
            end

            if badDataWarning
                warning('Bioinfo:geosoftread:BadGEOData',...
                    'Unable to read some lines of the file. Missing entries will be replaced with NaNs.');
            end
            outputStruct.Data = vals;
        catch theException %#ok<NASGU>
            warning('Bioinfo:geosoftread:IncompleteGEOFile','Problems reading the GEO data. The structure may be incomplete.');
        end

    case 'series'
        error('Bioinfo:geosoftread:SeriesNotSupported',...
            'Use GEOSERIESREAD to read in series data.');
    otherwise
        error('Bioinfo:geosoftread:UnknownGEOFormat',...
            'Unknown GEO format.\nGEOSOFTREAD supports Sample (GSM), Database (GDS) and Platform (GPL) records.');
end
