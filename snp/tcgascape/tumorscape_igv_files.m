function tumorscape_igv_files(input_parameter_file,run_dir,use_cap)
% make igv files for tumorscape

if ~exist('run_dir','var') || isempty(run_dir)
    run_dir = [pwd filesep];
end

% read tumorscape input parameter file (XML document)
TSP = read_tumorscape_params(input_parameter_file);

% cache file for cleaned D
cached_clean_master = [run_dir,'Dmaster.mat'];
D = load_D(cached_clean_master);

%% optionally cap the data
if exist('cap','var') && ~isempty(cap) && isnumeric(cap)
    % unravel parameter
    if isscalar(cap)
        maxcap = cap;
        mincap = -cap;
    elseif length(cap) == 2
        maxcap = cap(1);
        mincap = cap(2);
    else
        error('invalid length for cap parameter');
    end
    % cap the data
    if isa(D.dat,'SegArray')
        D.dat = cap_vals(D.dat,[maxcap mincap]);
    else
        D.dat(D.dat > maxcap) = maxcap;
        D.dat(D.dat < mincap) = mincap;
    end
end
%% save segfiles for IGV in local subdirectory
segfile_dir = [run_dir 'igvfiles/'];

if ~exist(segfile_dir,'dir');
    mkdir(segfile_dir);
    unix(['chmod 775 ' segfile_dir]);
end

write_igv_segs(D,'all_cancers',TSP.cancer_treegen,segfile_dir,TSP.cancer_namemap,40);
unix(['chmod 664 ' segfile_dir '/*']);

%% process segfiles for IGV

% matched IGV path and URL for served IGV directory
igv_dir = [TSP.igvfile_basedir filesep TSP.run_id];

%!!! temporary - use a different directory for testing
%! igv_dir = regexprep(igv_dir,'igvfiles','new_igvfiles');
%!!! end temporary

if ~exist(igv_dir,'dir');
    mkdir(igv_dir);
    unix(['chmod 775 ' igv_dir]);
end

%% create sample info file
xtypes = values(TSP.cancer_namemap,{D.sis.gcmtype});
SI = struct('Array',D.sdesc','Disease',xtypes);
infofilewrite([igv_dir '/sample_info.txt'],SI);

% gzip all the segfiles to the final IGV directory
d = dir(segfile_dir);
for f=1:length(d)
    fname = d(f).name;
    if length(fname) > 4 && strcmp(fname(end-3:end),'.seg')
        tname = [fname,'.gz'];
        unix(['gzip -c ',fullfile(segfile_dir,fname),'>',fullfile(igv_dir,tname)]); 
    end
end
unix(['chmod 664 ' igv_dir '/*']);

%% write version information
f = fopen([run_dir 'version.txt'],'w');
fprintf(f,'%s\n%s\n',tumorscape_version,gistic_version);
fclose(f);
unix(['chmod 664 ' run_dir '/version.txt']);

