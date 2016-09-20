function newfiles = get_new_filenames(oldfiles)
%NEWFILES = GET_NEW_FILENAMES(OLDFILES)
%  Finds a new filename for datastruct by incrementing number suffix on
%  filename with same root by 1.




if ~iscell(oldfiles)
    oldfiles = {oldfiles};
end


%oldfiles = regexprep(oldfiles,'.h5','');
oldfilesroot = regexp(oldfiles,[regexp_filesep '[^' regexp_filesep ']+\.'],'match');
oldfilesroot= regexprep([oldfilesroot{:}],regexp_filesep,'');

oldfilesroot = regexprep(oldfilesroot,'\.\d*\.','');
oldfilesroot = regexprep(oldfilesroot,'\.','');


prefix = regexp(oldfiles,[ '.+' regexp_filesep],'match'); 
fileprefix = [prefix{:}];

newfilenums = repmat({[]},1,length(oldfilesroot));

k=1;
for of = oldfilesroot

    try
        lst = ls([fileprefix{1} char(of) '*.h5']);  %try a listing of files with fileprefix
        lstcell = textscan(lst,'%s ')';
        lstcell = lstcell{:};
        lstcell = regexprep(lstcell,'.h5','');
        nums = regexp(lstcell,'\.\d*','match');
        nums(cellfun(@isempty,nums)) = {{'0'}};
        nums = [nums{:}];
        nums = regexp(nums,'\d*','match');

        nums = [nums{:}];
        nums = cellfun(@str2num,nums);
        nextnum = max(nums) + 1;
        nextnum = num2str(nextnum);
        newfilenums{k} = nextnum;
    catch
        newfilenums{k} = num2str(1);
    end
    k = k+1;
end
    


newfilesuffix = repmat({'.h5'},1,length(oldfiles));

newfiles = strcat(fileprefix,oldfilesroot,'.',newfilenums,...
    newfilesuffix);

