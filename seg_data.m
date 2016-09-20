%SEG_DATA - Loads submitted segment data for all centers

classdef seg_data
    properties
        basedir 
        cutoff_date
    end %properties
    methods
        function obj = seg_data()
            if ~ispc
                obj.basedir = '/xchip/tcga_scratch/DCC_downloads';
            else
                obj.basedir = '//thumper12/xchip_tcga_scratch/DCC_downloads';
            end
            obj.cutoff_date = 99999999;
        end
        function opendir = get_latest_open_dir(obj)
            d = dir ([obj.basedir '/20*-open']);
            b = sort(cellstr(char(d.name)));
            if length(b)<1
                throw(MException('seg_data:NoOpenDir','Could not find any open dirs'));
            end
            opendir_addon='';
            for i=length(b):-1:1
                if str2double(b{i}(1:8))<=obj.cutoff_date
                    opendir_addon = b{i};
                    break;
                end
            end
            if isempty(opendir_addon)
                throw(MException('seg_data:NoOpenDirMatching','Opendirs were found, but non before the cutoff date.'));
            end
            opendir = [obj.basedir '/' opendir_addon];
        end
        function restricteddir = get_latest_restricted_dir(obj)
            d = dir ([obj.basedir '/20*-restricted']);
            b = sort(cellstr(char(d.name)));
            if length(b)<1
                throw(MException('seg_data:NoRestrictedDir','Could not find any restricted dirs'));
            end
            restricteddir_addon='';
            for i=length(b):-1:1
                if str2double(b{i}(1:8))<=obj.cutoff_date
                    restricteddir_addon = b{i};
                    break;
                end
            end
            if isempty(restricteddir_addon)
                throw(MException('seg_data:NoRestrictedDirMatching','Restricteddirs were found, but non before the cutoff date.'));
            end
            restricteddir = [obj.basedir '/' restricteddir_addon];
        end

        function centerdir = get_center_dir(obj,center,batch)
            
            switch lower(center)
                case 'harvard'
                    centerdir_addon = ...
                        'tcga/tumor/ov/cgcc/hms.harvard.edu/hg-cgh-244a/cna';
                    centerdir = [obj.get_latest_open_dir,'/',centerdir_addon];
                case 'mskcc'
                    if 1 || batch~=9 %MSKCC has redone the batch 9 data on 1m
                        centerdir_addon = 'tcga/tumor/ov/cgcc/mskcc.org/cgh-1x1m_g4447a/cna';
                    else
                        centerdir_addon = 'tcga/tumor/ov/cgcc/mskcc.org/hg-cgh-244a/cna';
                    end
                    centerdir = [obj.get_latest_open_dir,'/',centerdir_addon];
                case 'stanford'
                    centerdir_addon = 'users/tcga4yeo/tumor/ov/cgcc/hudsonalpha.org/human1mduo/snp';
                    centerdir = [obj.get_latest_restricted_dir,'/',centerdir_addon];
                case 'broad'
                    centerdir_addon = 'users/tcga4yeo/tumor/ov/cgcc/broad.mit.edu/genome_wide_snp_6/snp';
                    centerdir = [obj.get_latest_restricted_dir,'/',centerdir_addon];
                    
                otherwise
                    ME = MException ('seg_data:UnrecognizedCenter',...
                        'Unrecognized center name.');
                    throw (ME);
            end

            
        end
        
        function platform = get_platform (obj,center,batch)
             switch (lower(center))
                case 'harvard'
                    platform = 'HG-CGH-244A';
                case 'mskcc'
                    if 1 || batch ~=9 % MSKCC has redone batch 9 on 1M
                        platform = 'CGH-1x1M_G4447A';
                    else
                        platform = 'HG-CGH-244A';
                    end
                case 'stanford'
                    platform = 'Human1MDuo';
                case 'broad'
                     platform = 'Genome_Wide_SNP_6';
                otherwise
                    ME = MException ('seg_data:UnrecognizedCenter',...
                        'Unrecognized center name.');
                    throw (ME);
            end
        end
        
        
        function batch_dir_name = get_batch_dir_name(obj,center, batch)
            
            switch lower(center)
                case 'harvard'
                    %TODO add in the alternative lookup used for GBM
                    harvard_ovarian_batch_lookup = [9 11 12 13 14 15];
                    if ~ismember(batch,harvard_ovarian_batch_lookup)
                        ME = MException ('seg_data:badHarvardbatchnum',...
                            'batch number not recognized in lookup table');
                        throw (ME);
                    end
                    center_batch = find(harvard_ovarian_batch_lookup==batch);
                case 'stanford'
                    %TODO add in the alternative lookup used for GBM
                    stanford_ovarian_batch_lookup = [9 14 15 17];
                    if ~ismember(batch,stanford_ovarian_batch_lookup)
                        ME = MException ('seg_data:badHarvardbatchnum',...
                            'batch number not recognized in lookup table');
                        throw (ME);
                    end
                    center_batch = find(stanford_ovarian_batch_lookup==batch);
                otherwise
                    center_batch = batch;
            end
            
            d = dir(obj.get_center_dir(center,batch));
            b = sort(cellstr(char(d.name)));
            %TODO technically, the .'s should be escaped in the regexp.
            %TODO the \d's should have a + qualifier after them
            switch (lower(center))
                case 'harvard'
                    re = ['^hms.harvard.edu_OV.' ...
                        obj.get_platform(center,batch) ...
                        '.' ...
                        num2str(center_batch) '.\d.\d$'];
                case 'mskcc'
%                     if batch ~=9
                        re = ['^mskcc.org_OV.' ...
                            obj.get_platform(center,batch) ...
                            '.' ...
                            num2str(center_batch) '.\d.\d$'];
%                     else
%                         re = ['^mskcc.org_OV.' ...
%                             obj.get_platform(center,batch) ...
%                             '.' ...
%                             num2str(center_batch) '.\d.\d$'];
%                     end
                case 'stanford'
                    re = ['^hudsonalpha.org_OV.' ...
                        obj.get_platform(center,batch) ...
                        '.' ...
                        num2str(center_batch) '.\d.\d$'];
                case 'broad'
                    re = ['^broad.mit.edu_OV.' ...
                        obj.get_platform(center,batch) ...
                        '.' ...
                        num2str(center_batch) '.\d.\d$'];
                otherwise
                    ME = MException ('seg_data:UnrecognizedCenter',...
                        'Unrecognized center name.');
                    throw (ME);
            end
            c = regexp(b,re);
            
            batch_dir_name_addon = '';
            for i=length(c):-1:1
                if isempty(c{i})
                    continue
                else
                    batch_dir_name_addon = b{i};
                    break
                end
            end
            if isempty(batch_dir_name_addon)
                ME = MException ('seg_data:badBatchnum',...
                    ['batch number' center_batch ...
                    ' not found in center''s data']);
                throw (ME);
            end
            
            batch_dir_name = [obj.get_center_dir(center,batch) '/' ...
                batch_dir_name_addon];
        end
        
        function batch_files = get_batch_files(obj,center, batch)
            switch lower(center)
                case 'harvard'
                    glob = '*_Segment.tsv';
                    d = dir([obj.get_batch_dir_name(center,batch) '/' glob]);
                    batch_filenames = {d(:).name};
                case 'mskcc'
                    if 1 || batch ~=9 % mskcc has redone batch9 on 1M
                        filename = ['mskcc.org_OV.CGH-1x1M_G4447A.' num2str(batch) '.CBS.txt'];
                    else
                        filename = ['mskcc.org_OV.HG-CGH-244A.' num2str(batch) '.CBS.txt'];
                    end
                    batch_filenames={filename};
                case 'stanford'
                    %TODO load non-tumor samples also...
                    dir_name = obj.get_batch_dir_name(center,batch);
                    dir_version_name = dir_name(end-4:end);
                    filename1 = ['hudsonalpha.org_OV.Human1MDuo.' dir_version_name '.seg.txt'];
                    filename2 = ['hudsonalpha.org_OV.Human1MDuo.' dir_version_name '.segnormal.txt'];
                    batch_filenames = {filename1,filename2};
                case 'broad'
                    glob = '*.seg.data.txt';
                    d = dir([obj.get_batch_dir_name(center,batch) '/' glob]);
                    batch_filenames = {'SampleInfo.txt' d(:).name};
                otherwise
                    ME = MException ('seg_data:UnrecognizedCenter',...
                        'Unrecognized center name.');
                    throw (ME);
            end
            %prepend directory name to each filename.
            batch_files = cellfun ...
                (@(x) horzcat(obj.get_batch_dir_name(center,batch),'/',x),...
                batch_filenames,...
                'UniformOutput',false);

        end
        
        function batch_data = get_batch_data(obj,center,batch)
            
            switch lower(center)
                case 'harvard'
                    batch_dir_name_tmp = ...
                        obj.get_batch_dir_name(center,batch);
                    batch_files = obj.get_batch_files(center,batch);
                    num_samples = size(batch_files,2);
                    if num_samples ==0 
                        ME = MException ('seg_data:NoSamples',...
                            ['No samples found in ' center ...
                            ' batch ' num2str(batch)]);
                        throw (ME);
                    end
                    %load all the samples
                    sample_values=cell(1,num_samples);
                    for i=1:num_samples
                        %print status
                        fprintf(1,'%d ',i);if ~mod(i,10) fprintf(1,'\n'); end
                        sample_values{i}=dataset(...
                            'file',batch_files{i},...
                            'format','%s%f%f%f%f','delimiter','\t',...
                            'treatasempty',{'NA','na'});
                    end
                    fprintf(1,'\n');
                    %add the sample name as a column, derived from filename
                    for i=1:num_samples
                        samplelen = size(sample_values{i},1);
                        sample_name= batch_files{i}(end-39:end-12);
                        sample_values{i}.Sample = ...
                            cellstr(repmat(sample_name,samplelen,1));
                    end
                        
                    %concatenate all in batch
                    batch_sample_values = vertcat(sample_values{:});
                    
                    %rearrange to standard column ordering
                    batch_sample_values = batch_sample_values(:,[6 1:5]);
                    
                    batch_data = batch_sample_values;
                case 'mskcc'
                    file = obj.get_batch_files(center,batch);
                    batch_data=dataset(...
                        'file',file{1},...
                        'format','%s%s%f%f%f%f%f','delimiter','\t',...
                        'treatasempty',{'NA','na'});
                    batch_data.num0x2Einformative=[];
                case 'stanford'
                    %TODO filter by batch number, rather than lumping all
                    %into one batch.
                    files = obj.get_batch_files(center,batch);
                    batch_data1=dataset(...
                        'file',files{1},...
                        'format','%s%s%f%f%f','delimiter','\t',...
                        'treatasempty',{'NA','na'});
                    batch_data2=dataset(...
                        'file',files{2},...
                        'format','%s%s%f%f%f','delimiter','\t',...
                        'treatasempty',{'NA','na'});
                    batch_data = [batch_data1;batch_data2];
                    %add empty probe number as column 5
                    batch_data.Probe_Number = NaN(size(batch_data,1),1);
                    batch_data = batch_data(:,[1:4 6 5]);
                case 'broad'
                    batch_dir_name_tmp = ...
                        obj.get_batch_dir_name(center,batch);
                    batch_files = obj.get_batch_files(center,batch);
                    num_files = size(batch_files,2);
                    num_samples = num_files-1;%first file is SampleInfo.txt file.
                    if num_samples < 1
                        ME = MException ('seg_data:NoSamples',...
                            ['No samples found in ' center ...
                            ' batch ' num2str(batch)]);
                        throw (ME);
                    end
                    %load sampleinfo file
                    file = batch_files{1};
                    sampleinfo = dataset(...
                        'file',file,...
                        'format','%s%s%s%s%s%s%s%s%s%s%s','delimiter','\t');
                    %load all the samples
                    sample_values=cell(1,num_samples);
                    for i=1:num_samples
                        file = batch_files{i+1};
                        %print status
                        fprintf(1,'%d ',i);if ~mod(i,10) fprintf(1,'\n'); end
                        sample_values{i}=dataset(...
                            'file',file,...
                            'format','%s%s%f%f%f%f','delimiter','\t',...
                        'treatasempty',{'NA','na'});
                    end
                    fprintf(1,'\n');
                    % Replace Broad sample name with TCGA sample name.
                    if batch < 14
                        for i=1:num_samples
                            broad_sample_name = sample_values{i}.ID(1);
                            % Look up TCGA name from sampleinfo
                            sampleind = ...
                                find (strcmp(sampleinfo.ARRAY,broad_sample_name));
                            if size(sampleind,1) < 1
                                throw(MException('seg_data:NameNotFound',...
                                    'Broad sample name not found in sample info list'));
                            end
                            if size(sampleind,1) > 1
                                throw(MException('seg_data:DuplicateName',...
                                    'Duplicate Broad sample name in sample info list'));
                            end
                            tcga_sample_name = sampleinfo.SAMPLE_NAME(sampleind);
                            % Replace first column with TCGA name
                            samplelen = size(sample_values{i},1);
                            sample_values{i}.ID = ...
                                cellstr(repmat(tcga_sample_name ,samplelen,1));
                        end
                    end
                    %concatenate all in batch
                    batch_data = vertcat(sample_values{:});
                    
                otherwise
                    ME = MException ('seg_data:UnrecognizedCenter',...
                        'Unrecognized center name.');
                    throw (ME);
            end
            
            %force uniform column naming convention
            batch_data = set (batch_data,'VarNames', {'Sample','Chromosome','Start','End',...
                'Num_Probes','Segment_Mean'});
            
            % Map x,y,m/mt, and xy to 23-26.
            inds = strcmpi(batch_data.Chromosome,'x');
            batch_data.Chromosome(inds)={'23'};
            inds = strcmpi(batch_data.Chromosome,'y');
            batch_data.Chromosome(inds)={'24'};
            inds = strcmpi(batch_data.Chromosome,'m');
            batch_data.Chromosome(inds)={'25'};
            inds = strcmpi(batch_data.Chromosome,'mt');
            batch_data.Chromosome(inds)={'25'};
            inds = strcmpi(batch_data.Chromosome,'xy');
            batch_data.Chromosome(inds)={'26'};

            %make the Chromosome column numeric
            batch_data.Chromosome = str2double(batch_data.Chromosome);
            
        end


        
    end %methods
    methods (Static)
       function probe_file = get_probe_file(center,batch)
           basedir = '/xchip/tcga_scratch/gsaksena/CancerGenomeAnalysisData/trunk/markerfiles/gistic_ovarian';
           switch lower(center)
               case 'harvard'
                   probe_file = 'harvard.probes.txt';
               case 'mskcc'
                   if 1 || batch~=9 %MSKCC has redone batch 9 on 1m.
                   probe_file = 'mskcc2.probes.txt';
                   else
                   probe_file = 'mskcc1.probes.txt';
                   end
               case 'stanford'
                   probe_file = 'hudsonalpha.probes.txt';
               case 'broad'
                   probe_file = 'broad.probes.txt';
               otherwise 
                   throw(MException ('seg_data:UnrecognizedCenter',...
                       'Unrecognized center name.'));
           end
           probe_file = [basedir '/' probe_file];
        end
        % These methods belong on their own separate object...
        function export_ds(ds_in, filename)
            export(ds_in,'file',filename);
        end
        function ds_out = filter_keep_tumors_only_ds(ds_in)
            samplenames = char(ds_in.Sample);
            samplecode = samplenames(:,14:15);
            inds = strmatch('01',samplecode,'exact');
            ds_out = ds_in(inds,:);
        end
        function ds_out = filter_keep_normals_only_ds(ds_in)
            samplenames = char(ds_in.Sample);
            samplecode = samplenames(:,14:15);
            inds = [strmatch('10',samplecode,'exact'); ...
                strmatch('11',samplecode,'exact')];
            ds_out = ds_in(inds,:);
        end
        function ds_out = filter_remove_excessively_segmented_samples_ds...
                (ds_in, iqr_threshold)
            if nargin<2
                iqr_threshold = 2;
            end
            if isempty(ds_in) %handle empty case
                ds_out=ds_in;
                return
            end
            grp_ds = grpstats(ds_in,'Sample','numel');
            grp_ds = sortrows(grp_ds,'GroupCount');
            gc = grp_ds.GroupCount;
            thresh = prctile(gc,75) + iqr_threshold*iqr(gc);
            inds_grp_keep = gc<thresh;
            samplenames_keep = grp_ds.Sample(inds_grp_keep);
            samples_keep = ismember(ds_in.Sample,samplenames_keep);
            ds_out=ds_in(samples_keep,:);
        end
        function ds_out = filter_remove_samples(ds_in, sample_blacklist)
            %sample_blacklist is a cell array of strings. If a string
            %matches the first part of a TCGA samplename, that sample will
            %be purged.
            samplenames = char(ds_in.Sample);
            dispose_inds = [];
            for i = 1:length(sample_blacklist)
                this_sample_name = sample_blacklist{i};
                sample_inds = strmatch(this_sample_name, samplenames);
                dispose_inds = [dispose_inds; sample_inds]; %#ok<AGROW>
            end
            ds_out = ds_in;
            ds_out(dispose_inds,:) = [];
        end
                
    end %methods(static)
end %classdef

    