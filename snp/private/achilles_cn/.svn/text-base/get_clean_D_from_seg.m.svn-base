function [CL SI] = get_clean_D_from_seg(cached_D_file,input,options)
% GET_CLEAN_D_FROM_SEG get/make pipelined copy number D structure from segments 

    if ~exist('input','var') && nargout > 1
        error('input parameter required.')
    end
    if ~isfield(input,'sampinfo_file') || isempty(input.sampinfo_file)
        error('valid sampinfo_file is a required input');
    else
        if ~exist(input.sampinfo_file,'file')
            error(['cannot find sample info file''',input.sampinfo_file,'''.'])
        else
            SI = read_sample_info_file(input.sampinfo_file);
        end
    end
    % see if data has already been cached
    if exist(cached_D_file,'file')
        
        %%% load cached preprocessed CN structure
        verbose('=> Reading saved CLE copy number data.',10);
        CL = load_D(cached_D_file);
        
    else
        %%% create cleaned CN structure and save it for next time
        if ~exist('input','var')
            error('input parameter required.')
        else
            if ~isfield(input,'markersfile') || isempty(input.markersfile)
                error('valid markers file is a required parameter');
            else
                if ~exist(input.markersfile,'file')
                    error(['cannot find markers file''',input.markersfile,'''.'])
                end
            end
            if ~isfield(input,'segfile') || isempty(input.segfile)
                error('valid segmented data file is a required parameter');
            else
                if ~exist(input.segfile,'file')
                    error(['cannot find segment file''',input.segfile,'''.'])
                end
            end
            % cnv file is optional
            input = impose_default_value(input,'cnv_file',[]);
        end
        % set defaults for optional paramters
        options = impose_default_value(options,'join_segment_size',10);
        options = impose_default_value(options,'use_segarray',true);
        options = impose_default_value(options,'remove_X',true);
        options = impose_default_value(options,'remove_nans',true);
        
        verbose('=> Constructing CL from CLE segment data.',10);
        CL = make_D_from_seg(input.segfile,input.markersfile,[],options.use_segarray);
        
        %% clean up segments
        verbose('=> Cleaning up data',10);
        
        % Remove CNVs
        if ~isempty(input.cnv_file) && exist(input.cnv_file,'file')
            verbose('Removing CNVs',20);
            [CL,nremoved] = remove_cnv(CL,input.cnv_file);
            fprintf(1,'Removed %d CNVs from raw data\r\n',nremoved);
        end
        % remove X,Y chromosome
        if options.remove_X
            verbose('Removing X chromosome...',20);
            CL=reorder_D_rows(CL,SegArray(CL.chrn)<=22); %! SegArray temp for memory optimization
        end
         
        % remove probes with a NaN value for any sample
        if options.remove_nans
            nans = any(isnan(CL.dat),2);
            if any(nans)
                verbose('removing %d probes w/NaN data',20,sum(nans));
                CL = reorder_D_rows(CL,~nans);
            end
        end
        
        % filter/smooth short segments
        CL.cbs = CL.dat;
        verbose('Merging small segments...',20);
        CL = smooth_cbs(CL,options.join_segment_size);
        CL.dat = CL.cbs;
        CL = rmfield_if_exists(CL,{'cbs','cbs_rl'});

        verbose('Median centering data...',20);
        CL.medians = median(CL.dat,1);
        CL.dat = CL.dat - repmat(CL.medians,size(CL.dat,1),1);

        % put cell-line information into SIS
        [~,matches] = match_string_sets_hash({SI.array},CL.sdesc);
        CL.sis = SI(matches);
        % remove some unused fields
        CL = rmfield_if_exists(CL,{'orig','origidx','gorigidx','chr','history'});
        % store information about how the data was created
        CL.create_info = struct;
        CL.create_info.date = datestr(now);
        CL.create_info.input = input;
        CL.create_info.options = options;
        
        % save the processed data
        save_D(cached_D_file,CL,'-v7.3');
    end
