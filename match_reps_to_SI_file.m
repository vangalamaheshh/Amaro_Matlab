function SI = match_reps_to_SI_file(SI,repnames,sifilename)
% MATCH_REPS_TO_SI_FILE - annotate a sample info file with duplicates found
% using FIND_DUPLICATE_SAMPLES.
%
%   SI = MATCH_REPS_TO_SI_FILE(SI,REPNAMES,SIFILENAME)
%
%       SI is the SI structure generated using read_sample_info_file
%
%       REPNAMES is the cell array of repnames generated using
%       FIND_DUPLICATE_SAMPLES
%
%       SIFILENAME is the name of the new sample info file.`   
%
%       A new column in the sifilename "duplicates_uniqids" lists the
%       uniqids of the samples that are duplicates.


siarrays = arrayfun(@(s) s.array,SI,'UniformOutput',0);
siuniqueids = arrayfun(@(s) s.uniqid,SI,'UniformOutput',0);
[SI.duplicates_uniqids] = deal('');
fldnms = fieldnames(SI);
SIcell = struct2cell(SI);
for k = 1:length(repnames)
    [m,i] = intersect(siarrays,repnames{k});
    alluniqids = unique(siuniqueids(i))';
    other_uniqids = cellfun(@(x,y) setdiff(x,y),repmat({alluniqids},length(i),1),siuniqueids(i),...
        'UniformOutput',0);
    
    semi = repmat({' ; '},1,length(alluniqids)-1);
    texttoadd = cellfun(@(x) reshape([x;semi],1,2*length(x)),other_uniqids,'UniformOutput',0);
    SIcell(strmatch('duplicates_uniqids',fldnms),i) = ...
        cellfun(@(x,y) cat(2,x,y),SIcell(strmatch('duplicates_uniqids',fldnms),i)',...
        cellfun(@(x) [x{:}],texttoadd,'UniformOutput',0),'UniformOutput',0);
end

SI = cell2struct(SIcell,fldnms,1);

if exist('sifilename','var')
    
    if exist(sifilename,'file')
        input('Filename already exists.  Overwrite file?')
        
        if 1==regexp(r,'^[yY]')
            infofilewrite(sifilename,SI);
        else
            warning('Not writing output.  Use ''infofilewrite'' to write SI'); %#ok
        end
        
    else
        infofilewrite(sifilename,SI)
    end
    
end

    