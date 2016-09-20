function [D,repnames,bad_tumor_names,arrays_for_core,celllines] = gistic_preprocessing(base_dir,...
    sample_info_file,array_list_file,Dscell,...
    snp_skip,show_hist,save_raw,norm_collapse_method,use_paired,...
    n_closest_n,perform_batch_correction,batch_effect_min_sz,batch_effect_bonf_pv_thresh,...
    batch_effect_abs_pv,output_dir_extension,hist_qc_normals,hist_qc_tumors,...
    conserve_memory,memlimit,use_all_cores,tan_ceil_val,scaleoper,startfromraw,genepattern,cnvfile,genomebuild)  %#ok<INUSL,INUSD>
%GISTIC_PREPROCESSING preprocess array data.  (Batch effect correction,
%normalization, and hist_qc.)
%
%   DOUT = GISTIC_PREPROCESSING(BASE_DIR, SAMPLE_INFO_FILE,
%   ARRAY_LIST_FILE, DFILESCELL, SNP_SKIP, SHOW_HIST, SAVE_RAW,
%   NORM_COLLAPSE_METHOD, USE_PAIRED, N_CLOSEST_N,
%   PERFORM_BATCH_CORRECTION, BATCH_EFFECT_MIN_SZ,
%   BATCH_EFFECT_BONF_PV_THRESH, BATCH_EFFECT_ABS_PV, OUTPUT_DIR_EXTENSION,
%   HIST_QC_NORMALS, HIST_QC_TUMORS, CONSERVE_MEMORY,MEMLIMIT,use_all_cores,tan_ceil_val,
%   scaleoper,startfromraw,genepattern,cnvfile,genomebuild)
%
%
%   Description of Inputs:
%
%   base_dir: (required) base directory for data processing; (string)
%   sample_info_file: (required) sample info file (in base dir) (string)
%   array_list_file: name of array list file (can be cell array of strings(?)), default empty (cell array of strings)
%   dscell: cell array listing filenames of D struct files (from
%       snp_to_D); pathnames are relative to base directory  OR cell array
%       of structs
%   snp_skip: decimation frequency on snp data, default 1 (int)
%   show_hist: show histogram?, 0 = no (default); 1 = yes (logical)
%   save_raw: save raw data?, 0 = no (default); 1 = yes (logical)
%   norm_collapse_method: 'mean', 'median', or 'tangent' (string)
%   use_paired: use_paired_norms? 0 = no (default); 1 = yes 
%   n_closest_n: how many norms for n_closest_n? default = 5 (int)
%   perform_batch_correction: 0 = no (default); 1 = yes (logical)
%   batch_effect_min_sz: minimum size batch to correct, default = 5 (int)
%   batch_effect_bonf_pv_thresh: pre-Bonferroni correction threshold,
%        default=0.05 (double) **NOTE: BONFERRONI CORRECTION NEVER
%        HAPPENS**
%   batch_effect_abs_pv: p value threshold for no Bonferroni correction,
%        default = 0.001 (double)
%   output_dir_extension: output dir extension (default = <yymmdd>)
%        (string)
%   hist_qc_normals: use hist_qc to remove normals with copy number changes?
%        (0 = no (default); 1 = yes) (logical)
%   hist_qc_tumors: use hist_qc to remove tumors without copy number changes?
%        (0 = no; 1 = yes (default))  (logical)
%   conserve_memory=1: uses shorts and uint8 to save memory
%                  =2: converts dstructs from structs to datastructs so that
%   data and affy_calls is saved as HDF5 files on disk.
%        (0 = no (default))
%   memlimit: when implemented in a subfunction, the maximum number of
%   bytes of data from a field in datastruct object to load into memory at
%   a time (should use conserve_memory=2). Default: 300000000 (300MB).
%   use_all_cores: use all cores for reference
%   tan_ceil_val: default empty
%   scaleoper: 'mean' by default.  Also accepts function handles.
%   startfromraw: 1 = load D.raw.mat rather than individual plates (Dscell should then have location of D.raw.mat file)
%   genepattern: 1 if running as executable, 0 if running in matlab (won't
%   datastructure files)
%   cnvfile (optional) the name of the cnvfile 
%   genomebuild: {hg16,hg17,hg18}
%  

%REMOVED:
%   save_after_tumor_qc: no = 0 (default); yes = 1.  (Save data before
%   removing bad tumors and samples not intended for gistic core?)
%   maxsamps, maxsnps
%
%   Change on 01 Feb 2008:  Selecting tangent normalization does not automatically cause histqc to be run on normals first. 
%  
% ---
% $Id$
% $Date: 2007-12-03 15:57:57 -0500 (Mon, 03 Dec 2007) $
% $LastChangedBy: rameen $
% $Rev$



%% Check Inputs

varlist1 = {'base_dir','sample_info_file','array_list_file','Dscell','snp_skip',...
    'show_hist','save_raw','norm_collapse_method','use_paired',...
    'n_closest_n','perform_batch_correction','batch_effect_min_sz','batch_effect_bonf_pv_thresh',...
    'batch_effect_abs_pv','output_dir_extension','hist_qc_normals','hist_qc_tumors','conserve_memory',...
    'memlimit','use_all_cores','tan_ceil_val','scaleoper','startfromraw','genepattern','cnvfile','genomebuild'};

defaults = {'ERR','ERR','ERR','{}','1','0','0','''median''','0','5','1','5','.05',...
    '.001','''output''' ,'0','1','0','300000000','1','[]','''mean''','0','0','[]','''hg18'''};

required = [1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

for idx = 1:length(varlist1)
    if required(idx) && (~exist(varlist1{idx},'var') || eval(['isempty(' varlist1{idx} ')']))
        error('Required input %s undefined.',varlist1{idx})
    elseif ~exist(varlist1{idx},'var') || eval(['isempty(' varlist1{idx} ')'])
        eval([varlist1{idx} '=' defaults{idx} ';'])
    end
end

if isempty(strmatch(filesep,fliplr(base_dir)))
    base_dir = [base_dir filesep];
end

%% Parse Normalization Method

if isstruct(norm_collapse_method)
    if strcmp(norm_collapse_method.method,'tangent')
        if isfield(norm_collapse_method,'norm_normals')
            tan_norm_normals = norm_collapse_method.norm_normals;
            if isfield(norm_collapse_method,'norm_atatime')
                atatime = norm_collapse_method.norm_atatime;
            else
                atatime = [];
            end
        else
            tan_norm_normals = 0;
        end
    end
    norm_collapse_method = norm_collapse_method.method;
end


%% Check existence and read/write permissions files and directories

%check permissions of base directory

if save_raw || conserve_memory==2
    if ~exist(base_dir,'dir')
        mkdir(base_dir);
    end

    [suc,perm] = fileattrib(base_dir);
    if suc ==0 || perm.UserWrite ~= 1 || perm.UserRead ~=1 || perm.UserExecute ~=1
        error('Insufficient access permission to %s',base_dir)
    end


    cd(base_dir)

    % Create Output directory
    if output_dir_extension(end) == filesep
        output_dir = [base_dir output_dir_extension];
    else
        output_dir = [base_dir output_dir_extension filesep];
    end

    mkdir(output_dir);
else

    if output_dir_extension(end) == filesep
        output_dir = [base_dir output_dir_extension];
    else
        output_dir = [base_dir output_dir_extension filesep];
    end

    mkdir(output_dir);
end

hdf5_dir = [base_dir 'hdf5' filesep];

if conserve_memory==2
  mkdir(hdf5_dir)

  verbose('Checking that HDF5 directory is empty',30)
  hdf5conts = dir(hdf5_dir);
  
  if length(hdf5conts)>2 % contains only . and ..
   error('HDF5 directory is not empty')
  end

end


% Create normalization output directory according to normalization type and
% check permissions

% 
% if hist_qc_normals
%     method_st = ['preprocess_' norm_collapse_method];
%     norm_dir=[ base_dir method_st '_' output_dir];
% 
%     if ~exist(norm_dir,'dir')
%         verbose(['Creating ' norm_dir],20);
%         mkdir(norm_dir);
%     else
%         warning('Normalization directory already exists');%#ok
%     end
% 
%     [suc,perm] = fileattrib(norm_dir);
% 
%     if perm.UserWrite ~= 1 || perm.UserExecute ~=1
%         error('Cannot write to directory ''%s''', norm_dir)
%     else
%         disp(['Write permissions OK for ',norm_dir])
%     end
% 
% end

% Make output_dir and check permissions
% if ~exist(output_dir,'dir')
%     verbose(['Creating ' output_dir],20);
%     mkdir(ouput_dir);
% else
%     warning('Output directory already exists');%#ok
% end
% 
% [suc,perm] = fileattrib(output_dir);
% 
% if perm.UserWrite ~= 1 || perm.UserExecute ~=1
%     error('Cannot write to directory ''%s''', output_dir)
% else
%     disp(['Write permissions OK for ',output_dir])
% end


% Check access permissions on sample_info_file, array_list_file, Dscell

[suc,perm] = fileattrib(sample_info_file);
if isempty(strmatch(filesep,sample_info_file))
    warning('Recommended to use full path name for SAMPLE_INFO_FILE input');  %#ok
end
if ~exist(sample_info_file,'file') || perm.UserRead ~=1
    error('Sample Info File does not exist or is unreadable.')
end


[suc,perm] = fileattrib(array_list_file);
if isempty(strmatch(filesep,array_list_file))
    warning('Recommended to use full path name for ARRAY_LIST_FILE input'); %#ok
end
if ~exist(array_list_file,'file') || perm.UserRead ~=1
    error('Array List File does not exist or is unreadable.')
end

if startfromraw
    hdf5dir = [base_dir 'hdf5'];
    mkdir(hdf5dir)
    D = load_D2(Dscell,hdf5dir);
else

    if ~iscell(Dscell)
        Dscell = {Dscell};
    end

    dclasses = cellfun(@class,Dscell,'UniformOutput',0);
    %if Dscell is list of filenames
    if ~isempty(strmatch('char',dclasses)) %if 1 is a filename, assume everything's a filename
        for k = 1:length(Dscell)
            if isempty(strmatch(filesep,Dscell{k}))
                warning('Recommended to use full path name for DFILESCELL input');%#ok
            end
            [suc,perm] = fileattrib(Dscell{k});
            if ~exist(Dscell{k},'file') || perm.UserRead ~=1
                error('Error with Dscell, file: %s;  File does not exist or is unreadable.',Dscell{k})
            end
        end
        use_ds_as_files = 1;
    else                    %assume structs or datastructs
        use_ds_as_files = 0;
    end

    verbose('File and directory permissions checked successfully.')


    %% Read sample info file
    SI = read_sample_info_file(sample_info_file);

    % add good field to SI as yes if it does not exist
    if ~isfield(SI,'good')
        [SI.good] = deal('EMPTY');
    end

    %% Read array list file

    AL = read_array_list_file(array_list_file);


    %% Match use array names to sample info

    [Mt,IDXalsi,IDXsial]=match_string_sets({AL.array},{SI.array});

    if length(IDXalsi)<length({AL.array}) || isempty(IDXalsi) % did not match all the arrays
       error('Did not find a match in the sample info file to all the arrays in the array list file');
    end

    
    SI = SI(IDXsial);
    AL = AL(IDXalsi);
    
    
    %check merge conditions OK
    verbose('Checking merge conditions',30);
    plats = {SI.platform};
    if isfield(AL,'inc_merge')
    inc_merge = cellfun(@str2num,{AL.inc_merge});
    uniqid = {SI.uniqid};
    [s,i,j] = unique(uniqid);


    multiplats  = find(histc(j,1:length(j))>1);
    illegal_merge = {};
    for k = multiplats
        idx = find(k==j);
        thisplats = plats(idx);
        mergeplats = thisplats(logical(inc_merge(idx)));

        if length(mergeplats) ~= length(unique(mergeplats))

            illegal_merge = [illegal_merge uniqid(idx(1))];%#ok
        end

    end
     end
%     if ~isempty(illegal_merge)
%         error('Illegal merge in unique ID: %s\n',illegal_merge{:});
%     end


    

    %% Load plates
    %----- Make sure array list file points to Dstructs. Array list file has
    %final say over who gets processed.
    verbose('Loading plates',20)

    allsamps = {AL.array};
    %
    if ~isempty(Dscell) && use_ds_as_files

        [AL,Dstructs]=load_Dstructs(Dscell,array_list_file,AL,conserve_memory,hdf5_dir,cnvfile,0);
        %make Dstructs cell array
        for k = 1:length(Dstructs)

            Dstructs{k} = setmetadata(Dstructs{k},'MemLimit',memlimit);
            allsamps = setdiff(allsamps,Dstructs{k}.sdesc);
     
        end

    elseif ~isempty(Dscell) && ~use_ds_as_files
        Dstructs = Dscell;

        if conserve_memory==1
            for k = 1:length(Dstructs)
                if isa(Dstructs{k},'struct')
                    if ~isa(Dstructs{k}.dat,'single')
                        Dstructs{k}.dat = single(Dstructs{k}.dat);
                    end
                    if ~isa(Dstructs{k}.affy_calls,'uint8')
                        Dstructs{k}.affy_calls(isnan(Dstructs{k}.affy_calls)) = 4;
                        Dstructs{k}.affy_calls = uint8(Dstructs{k}.affy_calls);
                    end

                      allsamps = setdiff(allsamps,Dstructs{k}.sdesc);

                else
                    error('must be a struct to use convserve_memory==1');
                end
                
                if exist('cnvfile','var') && ~isempty(cnvfile)
                    Dstructs{k} = remove_cnv(Dstructs{k},cnvfile,0);
                end
                Dstructs{k} = order_by_pos(Dstructs{k});
            end
        elseif conserve_memory==2
            for k = 1:length(Dstructs)
                if isa(Dstructs{k},'struct')
                    if ~isa(Dstructs{k}.dat,'single')
                        Dstructs{k}.dat = single(Dstructs{k}.dat);
                    end

                    if ~isa(Dstructs{k}.affy_calls,'uint8')
                        Dstructs{k}.affy_calls(isnan(Dstructs{k}.affy_calls)) = 4;
                        Dstructs{k}.affy_calls = uint8(Dstructs{k}.affy_calls);
                    end

                    Dstructs{k} = datastruct(Dstructs{k});
                    
                    if exist('cnvfile','var') && ~isempty(cnvfile)
                        Dstructs{k} = remove_cnv(Dstructs{k},cnvfile,0);
                    end
                    Dstructs{k} = order_by_pos(Dstructs{k});

                    dflds = intersect({'dat','affy_calls'},fieldnames(Dstructs{k}));
                    dfile = ['D' num2str(k)];
                    Dstructs{k} = convert_to_diskfield(Dstructs{k},dflds,strcat(hdf5_dir,dfile,'_',dflds,'.h5'),...
                        strcat(dflds,'_dataset'));

                elseif isa(Dstructs{k},'datastruct')

                    if exist('cnvfile','var') && ~isempty(cnvfile)
                        Dstructs{k} = remove_cnv(Dstructs{k},cnvfile,0);
                    end
                    Dstructs{k} = order_by_pos(Dstructs{k});
                    
                    dflds = setdiff({'dat','affy_calls'},diskfieldnames(Dstructs{k}));
                    dfile = ['D' num2str(k)];
                    Dstructs{k} = convert_to_diskfield(Dstructs{k},dflds,strcat(hdf5_dir,dfile,'_',dflds,'.h5'),...
                        strcat(dflds,'_dataset'));

                end
                allsamps = setdiff(allsamps,{Dstructs{k}.sdesc});
                Dstructs{k} = setmetadata(Dstructs{k},'MemLimit',memlimit);
            end
        else
            for k = 1:length(Dstructs)
                allsamps = setdiff(allsamps,Dstructs{k}.sdesc);
            end
            if exist('cnvfile','var') && ~isempty(cnvfile)
                Dstructs{k} = remove_cnv(Dstructs{k},cnvfile,0);
            end
            Dstructs{k} = order_by_pos(Dstructs{k});
        end

        clear Dscell;
    else
        %throw error: files for data structs must be given somewhere!!
        error('D structure file names must be given as inputs.')
    end

if ~isempty(allsamps)
    warning(['No data loaded for samples: ' repmat('%s ; ',1,length(allsamps))],allsamps{:});
end


    %% Build Data Structure (M):
    %------Trim Datasets (snp subsets and remove affy control)
    %------Add sample info and genome location

    for i=1:length(Dstructs)  %loop through data structures
        verbose(['Constructing D cell with Dstruct' num2str(i)],30)
        D = Dstructs{i};
        snp_subset=1:snp_skip:size(D.dat,1);
        %---- replace data with subset
        verbose('Reordering rows',30);
        D=reorder_D_rows(D,snp_subset);

        %---- remove AFFY control probesets
        %     affx_idx=grep('^AFFX',D.marker,1);
        %     D=reorder_D_rows(D,setdiff(1:length(D.marker),affx_idx));
        verbose('Removing ''orig'' field',30);
        D=rmfield_if_exists(D,'orig');
        %----- match samples to sample_info and add genome_info

        verbose('Adding sample info',30);
        D = add_sample_info(D,SI,'array');  %CHANGE THIS BACK to 'ARRAY'!!!

        verbose('Adding include information from AL file',30);
        D = add_array_includes(D,AL);
        verbose('Removing ''origidx'' field',30);
        D=rmfield_if_exists(D,'origidx');
        %-- sort according to genomic location - - move to before unite_Ds
%         verbose('Ordering D by pos',30);
%         D=order_by_pos(D);

        Dstructs{i} = D;
    end

    M = Dstructs;
    clear Dstructs

    %% Remove markers w/o genomic position and sort to genomic location
% 
%     for k = 1:length(M)
% 
%         w_pos=find(~isnan(M{k}.pos)); % --move to before unite_Ds
%         verbose(['Mraw(' num2str(k) '): Keeping ' num2str(length(w_pos)) ' markers with genomic position'],20);
%         M{k}=reorder_D_rows(M{k},w_pos);
% % 
% %         %% Sort according to genomic location  -- move to before unite_Ds
% %         verbose('Sorting according to genomic location',20);
% %         M{k}=order_by_pos(M{k});
%     end

    %% Sort Dstructs by platform and unite within platform

    M = sort_by_platform(M);
    for k = 1:length(M)
        markers = M{k}{1}.marker;
        for jj = 1:length(M{k})        
            markers = intersect(markers,M{k}{jj}.marker(~isnan(M{k}{jj}.pos))); %note that order is not preserved
        end
        for jj = 1:length(M{k})
            [cc,ia,ib] = intersect(markers,M{k}{jj}.marker);
            verbose('Sorting according to genomic location',20);
            M{k}{jj} = order_by_pos(reorder_D_rows(M{k}{jj},ib));
        end 
        
      
        M2{k} = unite_Ds(M{k}); %#ok

  
        for kk = 1:length(M{k})
            if strcmp('datastruct',class(M{k}{kk}))
            deleteDfiles(M{k}{kk});
            end
        end
    end

    if iscell(M) && length(M) == 1
        M2 = M2{1};
    end

    D = M2;
    clear M M2

    %% Save raw data (if needed)

    if save_raw
        verbose(['Saving raw data to: ' output_dir 'D.raw.mat'],20);

        D = save_D2([output_dir 'D.raw.mat'],D);  

    end

end


%     s = whos('M');
%     if s.bytes < 2000000000  %save without 'v7.3' tag is much faster; '-v7.3' only necessary for files over 9 GB
%         save([ output_dir 'M.D.mat'],'M');
%     else
%         save([ output_dir 'M.D.mat'],'M','-v7.3');
%     end
% % % % end
% % % % 
%
%  IMPLEMENT THIS WITH save_D2 function:
%
% % % % if iscell(M)
% % % %     for k = 1:length(M)
% % % %         if strmatch(class(M{k}),'datastruct','exact')
% % % %             M{k} = change_data_pointer(M{k},1);
% % % %         end
% % % %     end
% % % % elseif strmatch(class(M{k}),'datastruct','exact')
% % % %     M{k} = change_data_pointer(M{k},1);
% % % % end



% -v7.3 flag is slow!
% tic;save Mnoflag M;toc
% Elapsed time is 0.977114 seconds.
% tic;save Mwflag M -v7.3;toc
% Elapsed time is 537.302946 seconds.



%% Batch Effect 

D = preproc_log2trans(D,1,1);

D = preproc_scaledata(D,1,scaleoper);  %Remove while comparing to Gaddy's


if perform_batch_correction


    batch_effect_params=struct('method','batch_effect_correct','one_vs_all','yes',...
        'min_sz',batch_effect_min_sz,'bonf_pv_threshold',batch_effect_bonf_pv_thresh,...
        'absolute_pv',batch_effect_abs_pv);%,'cbediary_filename', [output_dir 'batcheffectdata.txt']);

%     save_D2([output_dir 'DbeforeBC.mat'],D); 

   % D = preproc_correctbatcheffect(D,batch_effect_params,output_dir);
    D = preproc_correctbatcheffect(D,batch_effect_params);
else
    verbose('Skipping Batch Effect', 20)
end


% 
% for k = 1:length(D)
%    w_pos=find(~isnan(D{k}.pos));
%     verbose(['D(' num2str(k) '): Keeping ' num2str(length(w_pos)) ' markers with genomic position'],20);
%     D{k}=reorder_D_rows(D{k},w_pos);
% end

%set collapse method to qc_collapse method if hist_qc_normals,
%otherwise should be norm_collapse_method



%-- We're merging before we do the histqc -- Q: is this what we want to do?
% A: Rameen says: Yes!

%% Normalization and quality control of normals.


if (hist_qc_normals) 
    if strcmp(norm_collapse_method,'tangent')
        qc_collapse_method = 'median';
    else
        qc_collapse_method = norm_collapse_method;
    end
    
    D = preproc_dividebynormals(D,n_closest_n,qc_collapse_method,use_paired,1);
    D = preproc_scaledata(D,1,scaleoper);  

    %% TO DO:  ADD SMARTER MERGE TO BUILD HDF5 file with each plate

        
    if length(D)>1     
        D = preproc_mergeplatforms(D);

        ismerged = 1;
    else
        ismerged = 0;
    end
    
    normdir = [output_dir 'normals_hqc'];
    D = preproc_histqcnormals(D,base_dir,normdir,1,show_hist);

    if isequal(norm_collapse_method,'tangent')
        verbose('Using tangent normalization.  Performing stage 2 histogram_qc on normals',10);
        
        if ismerged
            D = preproc_unmergeplatforms(D);
         
        end
        
        D = preproc_dividebynormals(D,n_closest_n,'tangent',use_paired,1);
        D = preproc_scaledata(D,1,scaleoper);  %Remove while comparing to Gaddy's

        if length(D)>1
            D = preproc_mergeplatforms(D);
        
        end
        
        D = preproc_histqcnormals(D,base_dir,normdir,2,show_hist);
    end

else
    
    if isequal(norm_collapse_method,'tangent')
       % D = preproc_fasttangent(D,use_all_cores,tan_ceil_val);
       if exist('tan_norm_normals','var') && tan_norm_normals
        D = preproc_slowtangent(D,atatime);
       else
        D = preproc_fasttangent(D);
       end
    else
        D = preproc_dividebynormals(D,n_closest_n,norm_collapse_method,use_paired,use_all_cores);
    end
    D = preproc_scaledata(D,1,scaleoper);  %Remove while comparing to Gaddy's

    

    if length(D)>1
     
        Dm = preproc_mergeplatforms(D);
        for k = 1:length(D)
            deleteDfiles(D{k});
        end
        clear D
       
        D = Dm;
        clear Dm;
  
    end
    
    verbose('Skipping histogram_qc on normals!',10)
   
    ploidy = get_sis(D,'ploidy');
    
    normals = strmatch({'2'},ploidy);
    v = zeros(size(D.dat,2),1);
    v(normals)=-1;
    %add a -1 to hist_qc in supdat to indicate no quality control on
    %normals
    D = add_D_sup(D,'hist_qc','histogram_quality_control',v','cols');

end



%% Quality Control on Tumors

if hist_qc_tumors
    verbose('Performing histogram_qc on tumors!',20);

    [one_peak, D] = histogram_qc(D,[],1,show_hist,1,0);

    normals = findsupdat(D,'N',1,'and');
    tumors = as_column(setdiff(1:size(D.dat,2),normals));
    bad_tumors = intersect(one_peak,tumors);

    if ndims(D.supdat ==3)
        v3 = D.supdat(find_supid(D,'hist_qc'),:,:);
        v3(:,bad_tumors,:)=3;
        D.supdat(find_supid(D,'hist_qc'),:,:) = v3;
    elseif ndims(D.supdat ==2)
        v3 = D.supdat(find_supid(D,'hist_qc'),:);
        v3(bad_tumors)=3;
        D.supdat(find_supid(D,'hist_qc'),:) = v3;
    end

    bad_tumor_names=D.sdesc(bad_tumors);%#ok
%     save([output_dir 'flag_bad_tumors.mat'],'bad_tumors','bad_tumor_names');
    
    good_tumors = as_column(setdiff(tumors,bad_tumors));
    
else
    verbose('Skipping histogram_qc on tumors!',20);
    normals = findsupdat(D,'N',1,'and');
    good_tumors = as_column(setdiff(1:size(D.dat,2),normals));
    bad_tumor_names = {'No quality control performed on tumors'};
end



% %add snp scores
if ~genepattern
    if ~exist('genomebuild','var') || isempty(genomebuild) || strcmp('hg18',genomebuild)
        load /xchip/gistic/variables/hg18/hg18_20080402.mat;
    elseif strcmp('hg17',genomebuild)
        load /xchip/gistic/variables/cyto_rg_hg17_ucsc20070227.mat;
    elseif strcmp('hg16',genomebuild)
        load /xchip/gistic/variables/hg16_20070112.mat;
    end
    D = add_D_snp_scores(D,cyto);
end

D = preproc_invlog2trans(D); 


% if ~save_after_tumor_qc
%     mfile = 'D_before_removal.mat';
%     save([output_dir mfile], 'D')
% end

%%  Remove bad tumors

if isempty(good_tumors)
    warning('No tumors passed quality control')
end
% 
% D = reorder_D_cols(D,good_tumors)  ;




%% Find Duplicates
repnames = [];
%%REMOVE THIS UNTIL MEM ISSUES FIXED%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [D,repids,repnames] = find_duplicate_samples(D,0.4,0);
% 
% remdups = [];
% 
% for kk = 1:length(repids)
%     verbose(strvcat(repnames{kk}),20);%#ok
%     keepsamps = repids{kk}(1);
%     remdups = [remdups setdiff(repids{kk},keepsamps)]; %#ok
% end
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% D = reorder_D_cols(D,setdiff(1:size(D.dat,2),remdups));
% % 
% save([output_dir 'removedduplicates'], 'remdups') ;

%% Select samples for core
% 
% D = check_array_includes(D,'core');
% 
% save([output_dir output_fname],'D')
% 
% if isa(D,'datastruct')
%     D = change_data_pointer(D,1);
% end
% 
% if ~isempty(remdups)
%     arrays_for_core = setdiff(D.sdesc(good_tumors),D.sdesc(remdups));
% else
     arrays_for_core = D.sdesc(good_tumors);
% end


%% Remove cell lines

if isfield(D.sis,'cellxeno')
    celllines = strmatch('cell line',get_sis(D,'cellxeno'));
else
    celllines = [];
end




verbose('Removing Cell Lines from Array List for Core',10);
arrays_for_core = intersect(arrays_for_core,D.sdesc(setdiff(1:length(D.sdesc),celllines)));

%% Save Data
if ~genepattern
    
    verbose(['Saving processing output to : ' output_dir ],20);
    
    
    
    % WRITE DUPLICATES
    
    %REMOVED UNTIL MEM ISSUES FIXED
%     dupsfile = [output_dir 'duplicates'];
% 
%     repnameschar = cellfun(@strvcat,repnames,'UniformOutput',0);
%     if ~isempty(repnames)
%         fid = fopen(dupsfile,'w');
%         for k = 1:length(repnames)
%           fprintf(fid,'Appear to be duplicates: \n');
%           for kk = 1:size(repnameschar{k},1)
%             fprintf(fid,'%s \n',repnameschar{k}(kk,:));
%           end
%           fprintf(fid,'\n');
%         end
%         fclose(fid);
%     else
%         fid = fopen(dupsfile,'w');
%         fprintf(fid,'No Duplicates Found');
%         fclose(fid);
%         verbose('No Duplicates Found',30);
%     end

%     % WRITE CELL LINES
%     celllinesfile = [output_dir 'cell_lines'];
%     if ~isempty(celllines)
%         fid = fopen(celllinesfile,'w');
%         fprintf(fid,'Cell Lines: \n')
%         for k = celllines
%             fprintf(fid,'%s \n',D.sdesc{k});
%         end
%         fclose(fid);
%     else
%         fid = fopen(celllinesfile,'w');
% 
%         fprintf(fid,'No Cell Lines Found');
% 
%         fclose(fid);
%         verbose('No Cell Lines Found',30);
%     end

    % WRITE BAD TUMORS
    badtumorsfile = [output_dir 'badtumors'];
    fid = fopen(badtumorsfile,'w');
    fprintf(fid,'Bad Tumors: \n');
    fprintf(fid,'%s \n',bad_tumor_names{:});
    fclose(fid);

    % WRITE OUTPUT ARRAY LIST
    if ~isempty(arrays_for_core)
        S = cell2struct(arrays_for_core,'array',1);
        arraylistout = [output_dir 'core_array_list_' datestr(now,'yymmdd') '.txt'];
        infofilewrite(arraylistout,S);
    end
    
    verbose(['Saving processed data to: ' output_dir 'D.mat'],20);

    D2 = save_D2([output_dir 'D.mat'],D);
    D = D2;
end

arrays_for_core = as_column(arrays_for_core);  % mokelly 081027.  Assumed as column downstream.  Will not be a column if there are no normal samples


