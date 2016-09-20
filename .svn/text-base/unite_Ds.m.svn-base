function D=unite_Ds(Ds,rc,check,fillmissingfields,initializewithD1)
%UNITE_DS unite the data structures contained in a cell array of data
%structures.
%
%   D = UNITE_DS(DS,RC,CHECK,FILLMISSINGFIELDS,INITIALIZEWITHD1) unites the data structures
%   in cell array DS.  RC toggles between column-wise unite (default or set
%   RC='col') in which data structures are concatenated along dim 2 (i.e.
%   markers are matched) or row-wise unite in which structures are
%   concatinated along dim 1 (set RC='row'). Warning: if used in 'row'
%   mode, sample names are NOT MATCHED, but an error can be thrown if CHECK
%   is 1.  IF FILLMISSINGFIELDS is set to 0, an error will be thrown when
%   the D or the sis is missing a field contained in the union across all
%   Ds.  If FILLMISSINGFIELDS is empty of 1, no error will be thrown, and
%   a default value will be used (either NaN or 0).
%
%   See also PREPROC_MERGEPLATFORMS.  
%
%   TO DO: 
%   1) modify to take union, rather than intersection, of fields  -done
%   2) modify to ask whether to check for inconsistencies in union of
%   fields, and if yes, to check for those inconsistencies
%   
%
%   Updated 17 Sept 07 to include .ref field
%   Updated 02 Oct 07 to make .sdesc field a row cell array rather than column
%   (did not fix case for if ischar(D.sdesc)). Jen Dobson (jdobson@broad)
%   Updated 18 Oct 07 with help documentation.
%   Updated 24 Oct 07 to remove fields that are not common to all Ds.
%   Updated 12 Dec 07 to write new D as datastruct object if any of input Ds are
%   datastructs.
%   Updated 27 Dec 07 -- Sorts the rows of the supdat for consistency across Ds.
%   (Before, supdat was merged without checking the supacc.)
%   Updated 18 Mar 08 -- Added support for gsup fields
%
%
%
%---
% $Id$
% $Date$
% $LastChangedBy$
% $Rev$


%  If any one of the input Ds is a datastruct, the output D will be a
%  datastruct


%% Get dimensions of matrix fields and check that data will combine

if exist('check','var') && strcmp('check',check)
    check = 1;
else
    check = 0;
end

dim1array = cellfun(@getsize,Ds,repmat({'dat'},1,length(Ds)),repmat({1},1,length(Ds)));
dim2array = cellfun(@getsize,Ds,repmat({'dat'},1,length(Ds)),repmat({2},1,length(Ds)));

if ~exist('initializewithD1','var') || isempty(initializewithD1)
    initializewithD1 = 0;
end

if ~exist('rc','var')
    rc='col';
end

if strncmpi(rc,'col',3)
    if dim1array == dim1array(1)*ones(1,length(Ds))
        dim1 = dim1array(1);
        dim2 = sum(dim2array);
    else
        error('Dimensions not consistent along dim 1')
    end
else
    if dim2array == dim2array(1)*ones(1,length(Ds))
        dim2 = dim2array(1);
        dim1 = sum(dim1array);
    else
        error('Dimensions not consistent along dim 2')
    end
end


%% Get fieldnames of the final D;
dflds = cellfun(@fieldnames,Ds,'UniformOutput',0);
dflds = cellfun(@strvcat,dflds,'UniformOutput',0);
dflds = cellstr(strvcat(dflds));  %#ok
dflds = unique(dflds);


%% Set class for D and initialize fields;

Ds_classes = cellfun(@class,Ds,'UniformOutput',0);

% if datastruct, get dskflds and hdf5dir; set dirtybit
isdatastruct = strmatch('datastruct',Ds_classes);
hdf5dir = [];
dskflds = {};
memlimit = Inf;
if ~isempty(isdatastruct)
    for k = isdatastruct'
        dskflds = [dskflds diskfieldnames(Ds{k})];  %#ok
        if isempty(hdf5dir)
            hdf5dir = regexp(get_datafile(Ds{k},'dat'),'/.+/','match');  %will use dir of first D as hdf5dir
        end
        thisml = getmetadata(Ds{k},'MemLimit');
        if ~isempty(thisml)
            memlimit = min(memlimit,thisml);
        end
    end
    dskflds = unique(dskflds);
    D = datastruct();
    %make diskfields
    D = add_diskfield(D,get_new_filenames(strcat(hdf5dir,dskflds,'.h5')),dskflds,repmat({[dim1 dim2]},1,length(dskflds)),cellfun(@get_D_field_props,dskflds,repmat({'class'},1,length(dskflds)),'UniformOutput',0));
    D = setmetadata(D,'DirtyBit',1);

    D = setmetadata(D,'MemLimit',memlimit);
    D = setmetadata(D,'ReadOnly',0);
end


memfields = setdiff(dflds,dskflds);  %dskflds is empty if not datastruct
memfields = setdiff(memfields,{'gorigidx','origidx','history'});  %don't merge these fields
%make memfields -- initialize using values given in get_D_field_props
for fl = memfields
    D.(char(fl)) = [];
    if initializewithD1

        if size(Ds{1},1) == dims1 && size(Ds{1},2) == 1
            cls = class(Ds{1}.(char(fl))(1));
            if ~strcmp(cls,{'cell','struct'})
                c = cast(0,cls);
                D.(char(fl)) = repmat(c,dim1,1);
            elseif strcmp(cls,'cell')  %#ok
                D.(char(fl)) = repmat({[]},dim1,1);
            elseif strcmp(cls,'struct')  %#ok
                D.(char(fl)) = repmat(struct(),dim1,1);
            end
        elseif size(Ds{1},1) == 1 && size(Ds{1},2) == dims2
            cls = class(Ds{1}.(char(fl))(1));
            if ~strcmp(cls,{'cell','struct'})
                c = cast(0,cls);
                D.(char(fl)) = repmat(c,1,dim2);
            elseif strcmp(cls,'cell')  %#ok
                D.(char(fl)) = repmat({[]},1,dim2);
            elseif strcmp(cls,'struct')  %#ok
                D.(char(fl)) = repmat(struct(),1,dim2);
            end
        elseif size(Ds{1},1) == dims1 && size(Ds{1},2) == dims2
            cls = class(Ds{1}.(char(fl))(1));
            if ~strcmp(cls,{'cell','struct'})
                c = cast(0,cls);
                D.(char(fl)) = repmat(c,dim1,dim2);
            elseif strcmp(cls,'cell')  %#ok
                D.(char(fl)) = repmat({[]},1,dim2);
            elseif strcmp(cls,'struct')  %#ok
                D.(char(fl)) = repmat(struct(),dim1,dim2);
            end
        elseif size(Ds{1},1) == dims1 && ischar(Ds{1}(1))  % assuming field is char array (like old sdesc)
            for k = 1:length(Ds)
                Ds{k}.(char(fl)) = cellstr(Ds{k}.(char(fl)));
            end
            D.(char(fl)) = {[]};


        else
            cls = class(Ds{1}.(char(fl))(1));
            if ~strcmp(cls,{'cell','struct'})
                c = cast(0,cls);
                D.(char(fl)) = c;
            elseif strcmp(cls,'cell')  %#ok
                D.(char(fl)) = {[]};
            elseif strcmp(cls,'struct')  %#ok
                D.(char(fl)) = struct();
            end
        end


    else
        if isequal(get_D_field_props(char(fl),'dims'),1)
            c = get_D_field_props(char(fl),'initialize');
            if ~isempty(c)
                D.(char(fl)) = repmat(c,dim1,1);
            else
                D.(char(fl)) = [];
            end
        elseif isequal(get_D_field_props(char(fl),'dims'),2)
            c = get_D_field_props(char(fl),'initialize');
            if ~isempty(c)
                D.(char(fl)) = repmat(c,1,dim2);

            else
                D.(char(fl)) = [];
            end
        elseif isequal(get_D_field_props(char(fl),'dims'),[1 2])
            c = get_D_field_props(char(fl),'initialize');
            if ~isempty(c)
                D.(char(fl)) = repmat(c,dim1,dim2);

            else
                D.(char(fl)) = [];
            end
        elseif isequal(get_D_field_props(char(fl),'dims'),0)
            c = get_D_field_props(char(fl),'initialize');
            if ~isempty(c)
                D.(char(fl)) = c;
            else
                D.(char(fl)) = [];
            end
        end
    end

end

%%%

%%

% 
% %first merge diskfields
% 
% for fl = dskflds
%     d1ct = 0;
%     d2ct = 0;
%     for k = 1:length(Ds)
%         addsz = getsize(Ds{k},char(fl));
%         D.(char(fl))(d1ct+1:addsz(1),d2ct+1:addsz(2)) = Ds{k}.(char(fl));
%         d1ct = d1ct + addsz(1);
%         d2ct = d2ct + addsz(2);
%     end
% end
% 
% 
% %remove HDF5 associated with the single Ds
% 
% for k = 1:length(Ds) 
%     deleteDfiles(Ds{k})
% end


%%

if ~exist('check','var') || isempty(check)
    check = 0;
end

if ~exist('fillmissingfields','var') || isempty(fillmissingfields)
    fillmissingfields = 1;
end

ufields = dflds'; 

if ~initializewithD1
    uprops = cellfun(@get_D_field_props,ufields,repmat({'unite'},1,length(ufields)),'UniformOutput',0);
    rowfields = ufields(strcmp('row',uprops));
    mtxfields = ufields(strcmp('matrix',uprops));
    colfields = ufields(strcmp('column',uprops));
    gsupfields = ufields(strcmp('gsup',uprops));    
    supfields = ufields(strcmp('sup',uprops));   
    unknownfields = ufields(cellfun(@isempty,uprops));  %#ok
else
    flds = fieldnames(Ds{1});
    rowfields = flds(cellfun(@(x) isequal(size(Ds{1}.(char(x))),[1 dim2]),flds));
    colfields = flds(cellfun(@(x) isequal(size(Ds{1}.(char(x))),[dim1 1]),flds));
    gsupfields = {'gsupdat','gsupacc','gsupdesc'};
    supfields = {'supdat','supacc','supdesc'};
    
    mtxfields = flds(cellfun(@(x) isequal(size(Ds{1}.(char(x))),[dim1 dim2]),flds));
    unknownfields = setdiff(flds,[rowfields' colfields' mtxfields']);   %#ok
end


if strncmpi(rc,'col',3)
    if check
        gf = @(X,f) X.(f);
    end

    for fl = colfields
        if check                %for each column field, make sure the field from each D is exactly the same
            b = cellfun(gf,Ds,repmat(fl,1,length(Ds)),'UniformOutput',0);
            pass = cellfun(@isequal,b,repmat({b{1}},1,size(b,2)));
            if sum(~pass)>0
                error('Not all fields equal for field ''%s''',char(fl))
            end
            D.(char(fl)) = b{1};
        else                    %if no check, just take the field from the first D that has it.
            for k = 1:length(Ds)
                D.(char(fl)) = Ds{k}.(char(fl));
                if ~isempty(D.(char(fl)))
                    break
                end
            end
        end
    end
    
    for fl = gsupfields
         if check                %for each column field, make sure the field from each D is exactly the same
            b = cellfun(gf,Ds,repmat(fl,1,length(Ds)),'UniformOutput',0);
            pass = cellfun(@isequal,b,repmat({b{1}},1,size(b,2)));
            if sum(~pass)>0
                error('Not all fields equal for field ''%s''',char(fl))
            end
            D.(char(fl)) = b{1};
        else                    %if no check, just take the field from the first D that has it.
            for k = 1:length(Ds)
                D.(char(fl)) = Ds{k}.(char(fl));
                if ~isempty(D.(char(fl)))
                    break
                end
            end
        end
    end 
    
    
    for fl = mtxfields
        d2ct = 0;
        for k = 1:length(Ds)
            addsz = [dim1array(k) dim2array(k)];
            if ~isfield(Ds{k},char(fl))
                if ~fillmissingfields
                    error('Missing field %s from plate number %d',char(fl),k)
                end
            else
                D.(char(fl))(1:addsz(1),d2ct+1:d2ct+addsz(2)) = Ds{k}.(char(fl));
            end
            d2ct = d2ct + addsz(2);
        end
    end

    for fl = rowfields
        
        d2ct = 0;
        for k = 1:length(Ds)
            addsz = dim2array(k);
            if ~isfield(Ds{k},char(fl))
                if ~fillmissingfields
                    error('Missing field %s from plate number %d',char(fl),k)
                end
            end
            try
                D.(char(fl))(d2ct+1:d2ct+addsz) = Ds{k}.(char(fl));
            catch   %if strings are stored as char array rather than cell array
                D.(char(fl))(d2ct+1:d2ct+addsz) = cellstr(Ds{k}.(char(fl)));
            end
            d2ct = d2ct + addsz;
         
        end
    end


    %now do sis, supdat,supdesc,supacc

    if isfield(D,'sis')
        
        gsf = @(X) fieldnames(X.sis);
        sf = cellfun(gsf,Ds,'UniformOutput',0);
        sf = cellfun(@strvcat,sf,'UniformOutput',0);
        sf = strvcat(sf);  %#ok
        sf = unique(cellstr(sf));

        siscell = repmat({[]},length(sf),dim2);
        D.sis = cell2struct(siscell,sf,1);
        
        d2ct = 0;
        
        for k = 1:length(Ds)
            missing = setdiff(sf,fieldnames(Ds{k}.sis));
            if ~isempty(missing) && ~fillmissingfields
                error('Missing field %s from sis of plate %d',missing,k)
            end
            if ~isempty(missing)
                for m = missing
                    [Ds{k}.sis.(char(m))] = deal('EMPTY');
                end
            end
            D.sis(1+d2ct:d2ct+length(Ds{k}.sis)) = orderfields(Ds{k}.sis);
            d2ct = d2ct + length(Ds{k}.sis);
        end
        
    end
    
    if isfield(D,'supdat')
        gsacc = @(X) X.supacc;
        allsupacc = cellfun(gsacc,Ds,'UniformOutput',0);
        allsupacc = strvcat(allsupacc);  %#ok
        allsupacc = cellstr(allsupacc);
        allsupacc = unique(allsupacc);
        D.supacc = strvcat(allsupacc);  %#ok
        D.supdesc = D.supacc;
        D.supdat = zeros(size(allsupacc,1),dim2);
        
        d2ct = 0;
        
        for k = 1:length(Ds)
            dssupacc = cellstr(Ds{k}.supacc);
            [tf,srt] = ismember(allsupacc,dssupacc);
            b = nan(length(tf),size(Ds{k}.supdat,2));
            b(tf,:) = Ds{k}.supdat(srt(tf),:);
            D.supdat(:,d2ct+1:d2ct+size(Ds{k}.supdat,2)) = b;
            d2ct = d2ct + size(Ds{k}.supdat,2);

        end
        
    end
    

    
    
elseif strncmpi(rc, 'row',3)
    
    for fl = rowfields
        if check            %if check, make sure all row fields are exactly the same (columns must be ordered before calling this function)
            gf = @(X,f) X.(f);
            b = cellfun(gf,Ds,repmat(fl,1,length(Ds)),'UniformOutput',0);
            pass = cellfun(@isequal,b,repmat({b{1}},1,size(b,2)));
            if sum(~pass)>0
                error('Not all fields equal for field ''%s''',char(fl))
            end
            D.(char(fl)) = b{1};
        else                %if no check, just set the field equal to the field of the first D that has it.
            for k = 1:length(Ds)
                D.(char(fl)) = Ds{k}.(char(fl));
                if ~isempty(D.(char(fl)))
                    break
                end
            end
        end
    end

    for fl = mtxfields
        d1ct = 0;
        for k = 1:length(Ds)
            addsz = [dim1array(k) dim2array(k)];
            if ~isfield(Ds{k},char(fl))
                if ~fillmissingfields
                    error('Missing field %s from plate %d',char(fl),k)
                end
            else     
                D.(char(fl))(d1ct+1:d1ct+addsz(1),1:addsz(2)) = Ds{k}.(char(fl));
            end
            d1ct = d1ct + addsz(1);
        end
    end

    for fl = colfields
        d1ct = 0;
        for k = 1:length(Ds)
            addsz = dim1array(k);
            if ~isfield(Ds{k},char(fl))
                if ~fillmissingfields
                    error('Missing field %s from plate %d',char(fl),k);
                end
            else
                try
                    D.(char(fl))((d1ct+1):(d1ct+addsz)) = Ds{k}.(char(fl));
                catch   %if strings are stored as char array rather than cell array
                    D.(char(fl))((d1ct+1):(d1ct+addsz)) = cellstr(Ds{k}.(char(fl)));
                end


            end

            d1ct = d1ct + addsz;
        end
    end
           
    for fl = supfields

         if check                %for each column field, make sure the field from each D is exactly the same
            b = cellfun(gf,Ds,repmat(fl,1,length(Ds)),'UniformOutput',0);
            pass = cellfun(@isequal,b,repmat({b{1}},1,size(b,2)));
            if sum(~pass)>0
                error('Not all fields equal for field ''%s''',char(fl))
            end
            D.(char(fl)) = b{1};
        else                    %if no check, just take the field from the first D that has it.
            for k = 1:length(Ds)
                D.(char(fl)) = Ds{k}.(char(fl));
                if ~isempty(D.(char(fl)))
                    break
                end
            end
        end
    end 


    %Merge SIS
    if isfield(D,'sis') 
                D.sis = Ds{1}.sis;
                if check
                checknames = fieldnames(D.sis);
                end
        for k = 2:length(Ds)
            if check
                checknames = intersect(checknames,fieldnames(Ds{k}.sis));
                for fl = checknames'
                    if ~isequal(get_sis(D,char(fl)),get_sis(Ds{k},char(fl)))
                        error('SIS mismatch at plate %d',k)
                    end
                end
            end
            
            missing = setdiff(fieldnames(Ds{k}.sis),fieldnames(D.sis))';
            if ~isempty(missing)
                if ~fillmissingfields
                    error('Missing sis fields %s in plates before %d',missing,k)
                else
                    [D.sis.(char(missing))] = Ds{k}.sis.(char(missing));
                end
            end
        end
    end
    
    %Merge supdat,supdesc,supacc
    if isfield(D,'gsupdat')
        gsacc = @(X) X.gsupacc;
        allgsupacc = cellfun(gsacc,Ds,'UniformOutput',0);
        allgsupacc = strvcat(allgsupacc);  %#ok
        allgsupacc = cellstr(allgsupacc);
        allgsupacc = unique(allgsupacc);
        D.gsupacc = strvcat(allgsupacc);   %#ok
        D.gsupdat = zeros(size(allgsupacc,1),dim1);
        D.gsupdesc = D.gsupacc;
        
        d1ct = 0;   %#ok
        
        for k = 1:length(Ds)
            dsgsupacc = cellstr(Ds{k}.gsupacc);
            [tf,srt] = ismember(allgsupacc,dsgsupacc);
            b = nan(length(tf),size(Ds{k}.gsupdat,2));
            b(tf,:) = Ds{k}.gsupdat(srt(tf),:);

            D.gsupdat(tf,(d1ct+1):(d1ct+size(Ds{k}.gsupdat,2))) = b;
            d1ct = d1ct + size(Ds{k}.gsupdat,2);            
            
        end
        
    end
   
end

D = setmetadata(D,'DirtyBit',0);
       
