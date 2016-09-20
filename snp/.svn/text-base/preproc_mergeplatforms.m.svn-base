 function Cout = preproc_mergeplatforms(Cin)
%
%       COUT = PREPROC_MERGEPLATFORMS(Cin)
%
%         
%
%       Array List File Include Behaviour: 
%                   -1 : removed from D
%                    0 : forced NaN merge
%                    1 : merge as normal, with NaNs if necessary
%                    2 : merge as normals, with NaNs if necessary
%
%
%       Revisions:
%
%               3 Dec 07 -- added Cout has same class as Cin{1} (datastruct
%               or struct)
%
%
%to include: boolean to ask if want to throw out unmatched arrays
%also: if matching of arrays is not one-to-one (i.e. 2 xba for 1 hind) need
%to decide what to do.
%
%TO DO: loop the writing of adddat to stay within memory limits


%% Handle array list include information and get arrays to merge

procname = 'merge';

 
numCs = length(Cin);

for k = 1:numCs   %get list of uniqids that merge for each platform
    %
    %     mergenames{k} = get_sis(Cin{k},'uniqid');
    %
    %     wontmerge{k} = find(Cin{k}.supdat(strmatch('inc_merge',Cin{k}.supacc),:)==0);
    %     mergenames{k}(wontmerge{k}) = strcat(mergenames{k}(wontmerge{k}),repmat({['_' num2str(k) '*']},1,length(wontmerge{k})));
    %     mergenames{k} = mergenames{k}(find(~cellfun(@isempty,mergenames{k})));
    %     mergenames{k} = mergenames{k}(find(~strcmp(mergenames{k},'EMPTY')));

  
   if isfield(Cin{k},'supdat')
Cin{k} = check_array_includes(Cin{k},'merge');
 end
  
    mergenames{k} = get_sis(Cin{k},'uniqid');
    trnc = @(x) strtok(x,'###');
    mergenames{k} = cellfun(trnc,mergenames{k},'UniformOutput',0);
    if isfield(Cin{k},'supdat')
    wontmerge{k} = find(Cin{k}.supdat(strmatch('inc_merge',Cin{k}.supacc),:)==0);
    else
      wontmerge{k} = [];
    end
    
    mergenames{k} = mergenames{k}(setdiff(1:length(mergenames{k}), ...
                                          wontmerge{k}));
    
      
end


commonnames = mergenames{1};   %make a list of the union of all mergenames
for k = 2:numCs
    commonnames = union(commonnames, mergenames{k});
end


%% Get the samples in the same order.  Fill unmatched samples with NaN


lencomm = length(commonnames);

mkuarray = {};  %indices to give to kth struct to get in commonnames order

needsNaN = zeros(numCs,lencomm);  %boolean matrix to indicate which samples are missing and need NaNs

%Find the fields to merge
usefields = [];
for k = 1:numCs
    usefields = [usefields; fieldnames(Cin{1})];  %#ok
end
usefields = unique(usefields);
diskfields = usefields(strcmp('disk',cellfun(@get_D_field_props,usefields',repmat({'storage'},1,length(usefields)),'UniformOutput',0)))';

for k = 1:numCs
    
    [uniqnames,i,j] = unique(mergenames{k});  %check to make sure all names in platform are uniq        
    if length(uniqnames)<length(mergenames{k})  %throw error if uniqid's are not uniq across platform
        nonuniq = mergenames{k}(setdiff([1:length(mergenames{k})],i));
        error(['Unique ID %s of array ' num2str(k) ' is not unique!\n'],nonuniq{:});
    end
    
    [Match,mku,muu] = match_string_sets(mergenames{k},commonnames);

    if length(mku) < lencomm
        missing_from_k = setdiff([1:lencomm],muu);
        mergenames{k} = [mergenames{k} commonnames(missing_from_k)];
        needsNaN(k,[length(mku)+1:lencomm]) = 1;
        
        %get new mku
        [Match,mku,muu] = match_string_sets(mergenames{k},commonnames);

        %Fill in new names
        Cin{k}.sdesc = mergenames{k};
        
        %put in NaNs

        %Extend the size of the diskfields
        datsize = getsize(Cin{k},'dat');
        if strcmp('datastruct',class(Cin{k})) %&& isdiskfield(Cin{k},'dat')
            for fl = diskfields
                if ~exist('hdf5dir','var')
                    hdf5dir = regexp(get_datafile(Cin{k},char(fl)),[regexp_filesep '.+' regexp_filesep],'match');
                end
                
                Cin{k} = add_diskfield(Cin{k},[char(hdf5dir) 'Cm' num2str(k) '_ext' char(fl) '.h5'],['ext' char(fl)],...
                    getsize(Cin{k},char(fl)) + [0 length(missing_from_k)],get_D_field_props(char(fl),'class'));  %make an extended field
                Cin{k}.(['ext' char(fl)])(:,1:getsize(Cin{k},char(fl),2)) = Cin{k}.(char(fl));
                Cin{k} = changefieldname(Cin{k},char(fl),['old' char(fl)]);  %keep olddat around so that file is deleted on cleanup
                Cin{k} = changefieldname(Cin{k},['ext' char(fl)],char(fl));  %make replace the field with the extended field
                if isfloat(Cin{k}.(char(fl)))  %fill missing samples with NaNs                  
                    Cin{k}.(char(fl))(:,(getsize(Cin{k},['old' char(fl)],2)+1):(getsize(Cin{k},['old' char(fl)],2)+length(missing_from_k))) = NaN;
                elseif isinteger(Cin{k}.(char(fl)))
                    Cin{k}.(char(fl))(:,(getsize(Cin{k},['old' char(fl)],2)+1):(getsize(Cin{k},['old' char(fl)],2)+length(missing_from_k))) = intmax(class(Cin{k}.(char(fl))));  %if uint8, use 255 for NaN
                else
                    error('Diskfield data is not int or float')
                end
            end
        end


        if isfield(Cin{k},'sis')  %EMPTY FILL SIS
            sisfield = Cin{k}.sis;
            nansis = sisfield(1);
            for fl = fieldnames(nansis)'
                nansis.(char(fl)) = 'NaN';
            end
            sisfield(end+1:end+length(missing_from_k)) = nansis;    
            Cin{k}.sis = sisfield;
        end

        if isfield(Cin{k},'supdat') %NaN fill supdat
            Cin{k}.supdat(:,end+1:end+length(missing_from_k)) = NaN;
        end
        
    end


    mkuarray{k} = mku;  %the indices mapping Cin{k} to the new order given by commonnames 



end



%% Perform Merge
% M is merged Data struct
% Initialize M and set class equal to class of Cin{1}    
M = [];

if ~isempty(strmatch('datastruct',cellfun(@class,Cin,'UniformOutput',0)))
    M = datastruct(M);
    M = setmetadata(M,'DirtyBit',1);
    mlimt = getmetadata(Cin{1},'MemLimit');
    if isempty(mlimt)
        mlimt = 300000000;
    end
    
    M = setmetadata(M,'MemLimit',mlimt);
    for k = 1:numCs
        M = setmetadata(M,'MemLimit',min(getmetadata(Cin{k},'MemLimit'),getmetadata(M,'MemLimit')));
    end
    M = setmetadata(M,'ReadOnly',0);
end

% Initialize the fields of M
for field = usefields'
    M.(char(field)) = [];
end


%Rename the merged samples
M.sdesc = cellstr(commonnames)';

%Get dimensions of M fields
datdim2 = length(M.sdesc);
lenfun = @(x) size(x.dat,1);
nmarkers = cellfun(lenfun,Cin);
datdim1 = sum(nmarkers);

%initialize disk fields ('dat','affy_calls') if datastruct
if isa(M,'datastruct') 
    
    for fl = diskfields

        M = add_diskfield(M,[char(hdf5dir),char(fl),'_merged','.h5'],char(fl),[datdim1 datdim2],get_D_field_props(char(fl),'class'));
    end
   
end
    


%Get the genotype supacc to merge
gsupaccs = [];
gsupdescs = [];
for k = 1:numCs
    if isfield(Cin{k},'gsupacc')  %first, get rid of 'PLAT' info, since it might change
        Cin{k} = reorder_D_supdat(Cin{k},'row',setdiff([1:size(Cin{k}.gsupdat,1)],strmatch('PLAT',Cin{k}.gsupdat)));
        if isfield(Cin{k},'gsupacc') && ~isempty(Cin{k}.gsupacc)
            gsupaccs = [gsupaccs cellstr(Cin{k}.gsupacc)'];
            gsupdescs = [gsupdescs cellstr(Cin{k}.gsupdesc)'];
        end
    end
end

[ugsupaccs,ugidx,dum] = unique(gsupaccs);

ugsupdescs = gsupdescs(ugidx);

M.gsupdat = zeros(length(ugsupdescs)+1,datdim1); %add one for platform field
%M.dat = zeros(datdim1,datdim2);  %we should be preallocating this space

 
    
platformkey = [];

%%
% loop through the Cin




% merge disk fields first
%TO DO: loop the writing of adddat to stay within memory limits
for field = diskfields
    lastrow = 0;
    for k = 1:numCs
        sortkuidx = mkuarray{k};
        adddat = Cin{k}.(char(field));
        M.(char(field))(lastrow+1:lastrow+size(adddat,1),:) = adddat(:,sortkuidx);

        lastrow = lastrow + size(adddat,1);
    end
end


for k = 1:numCs

    for field = intersect(usefields,{'marker','pos','chrn','chr'})
        old = M.(char(field));
        add = Cin{k}.(char(field));
        if size(add,1) < size(add,2)  %make sure column data
            add = add';
        end
        M.(char(field)) =[old; add];
    end

    % put original platform number info in the gsup
    if isfield(Cin{k},'gsupacc') && ~isempty(Cin{k}.gsupacc)
        [mgs,igs,jgs] = match_string_sets(strtok(cellstr(Cin{k}.gsupacc),':'),ugsupaccs);
        missinggsups = setdiff([1:length(ugsupaccs)],jgs);
        if ~isempty(missinggsups)
            warning('Platform %d: Missing gsupdat for fields: %s \n\tAdding NaNs.', k,char(ugsupaccs(missinggsups))')
        end
        jgs = [jgs missinggsups];
        %add gsup for missing sups, fill with NaNs, throw warnings; append
        %missing sups to jgs
        tmpgsupdat = [Cin{k}.gsupdat; NaN(length(missinggsups),size(Cin{k}.gsupdat,2))];
    else
        tmpgsupdat = NaN(length(ugsupdescs),nmarkers(k));
    end

    M.gsupdat(:,sum(nmarkers([1:k]))-nmarkers(k)+1:sum(nmarkers([1:k]))) = [tmpgsupdat ; k.*ones(1,nmarkers(k))];

    if isfield(Cin{k},'sis')
        sisfield = Cin{k}.sis;
        thisplat = sisfield(1).platform;
    else
        thisplat = ['unknown' num2str(k)];
    end


    if isfield(M,'used_normals')
        un =Cin{k}.used_normals;
        M.used_normals.(['p' thisplat]) = un(sortkuidx);
    end

    platformkey = [platformkey thisplat ' -- ' num2str(k) ' || '];



end

M.gsupacc = strvcat([gsupaccs 'PLAT']);
M.gsupdesc = strvcat([gsupdescs ['Platform:' platformkey]]);

%%

if isfield(M,'sis')
    M.sis = preproc_mergesis(Cin,mkuarray);
end

if isfield(M,'supdat')
    [M.supacc,M.supdesc,M.supdat,M.mergsupdat] = preproc_mergesup(Cin,mkuarray);   %should we just store the supdat etc and then merge once we're out of the loop?
end



%straight copy over any fields that were not merged
for fl = fieldnames(M)'
    if isempty(M.(char(fl)))
        M.(char(fl)) = Cin{1}.(char(fl));
    end
end

if isfield(M,'origidx')
    M = rmfield(M,'origidx');
end
if isfield(M,'gorigidx')
    M = rmfield(M,'gorigidx');
end

if isfield(M,'origidx')
    M = rmfield(M,'origidx');
end

if isfield(M,'ref')
    M = rmfield(M,'ref');
end

for k = 1:length(Cin)
    deleteDfiles(Cin{k});
end




%% Sort according to genomic location

M=reorder_D_rows(M,find(~isnan(M.pos)));

M=order_by_pos(M);

M = setmetadata(M,'DirtyBit',0);
Cout = M;
