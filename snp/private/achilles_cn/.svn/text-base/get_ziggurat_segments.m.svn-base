function [Qs D] = get_ziggurat_segments(cached_Qs_file,D,cyto,options)
    
    options = impose_default_value(options, 'max_segs_per_sample', 5000);
%   options = impose_default_value(options, 'sample_selection',1:size(D.dat,2));
    options = impose_default_value(options, 'cn_cap',1.5);
    
    ziggs = struct ; % ziggurat deconstruction options
    ziggs.max_segs_per_sample = options.max_segs_per_sample;
    if ~isfield(D,'islog')
        D.islog = true;
    end
    if isfield(D,'Qs')
        Qs = D.Qs;
    else
        % if the ziggurat deconstruction has been cached,
        % set the option to just read the data from a file
        if exist(cached_Qs_file,'file')
            verbose('loading ziggurat results from file ''%s''.',20,cached_Qs_file);
            Qs = load_Qs(cached_Qs_file);
            D.Qs = Qs;
        else
            verbose('performing ziggurat analysis of CN data',10);
            [D,Qs] = perform_ziggurat_deconstruction(D,cyto,ziggs,options.cn_cap);
            if ~D.islog
                D.dat = log2(D.dat + 2) - 1;
                D.islog = true;
            end
        end
        if isfield(options,'D_file')
            verbose('saving D.Qs in ''%s''.',20,options.D_file); 
            save_D(options.D_file,D,'-v7.3');
        end
    end
