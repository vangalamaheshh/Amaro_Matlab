function gp_gistic_preprocessing(varargin)
%
% gp_gistic_preprocessing 
%-o: (string) output directory (REQUIRED) 
%-si sample_info_file (REQUIRED) 
%-al array_list_file (REQUIRED)
%-od output_dir_ext 
%-of output_file_name  (REMOVE) 
%-ns norm_selection
%-nc norm_collapse_type 
%-ncn n_closest_n 
%-up use_paired_norm 
%-bc correct_batch_effect 
%-bems batch_effect_min_sz 
%-bebpv batch_effect_bonf_pv 
%-beapv batch_effect_absolute_pv 
%-ss snp_skip 
%-svr save_raw 
%-hn hist_qc_normals 
%-ht hist_qc_tumors 
%-sh show_hist 
%-cm conserve_memory
%-sm save_mat  = save D as matlab structure?
%-mfn maxfornorm = maximum number of arrays to use as Y in dist function
%-v verbose level
%
%      The options passed to GP_SNP_TO_D are:
%-i: (string) input snp file name (REQUIRED) 
%-a: (int) allele data (0 if no allele data; 1 if allele data as [intensity 1; call; intensity 2], 2 if allele data given as [intensity 1; intensity 2; call] -nl: (int) number of lines in file 
%-ec: (int) number of extra columns in file between chromosome position column and intensity data
%-h: number of header lines
%-sz: number of lines of data to read at each iteration of TEXTSCAN
%
%
%
% input: a directory of .D.mat files (for the various plates) 
%        
%        a master sample_info file
%        a list of tumor and normal samples to use
%       
%  COMPILE INSTRUCTIONS:
%
%       1) run 'mbuild -setup'.  Verify that the mbuildopts.sh that gets
%       installed in ~/.matlab/R2007b references the gcc in /util/gcc-4.1.1
%        (in instances where CC is set to 'gcc', make sure that the gcc is
%        '/util/gcc-4.1.1/bin/gcc'.)
%
%       2) go to the directory where you want to compile and type:
%           mcc -v -m -w enable gp_gistic_preprocessing
%
%       3) After compilation, you are ready to execute.  Set the
%       LD_LIBRARY_PATH environment variable by typing:
%           'setenv LD_LIBRARY_PATH
%           /broad/tools/apps/matlab75/sys/os/glnxa64:/broad/tools/apps/ma
%           tlab75/bin/glnxa64:/util/gcc-4.1.1/lib:$LD_LIBRARY_PATH'
%
%       4) execute the compiled binary
% ---
% $Id$
% $Date: 2007-12-03 15:57:57 -0500 (Mon, 03 Dec 2007) $
% $LastChangedBy: rameen $
% $Rev$
  
verbose([datestr(now) ':Running GP_GISTIC_PREPROCESSING'])

%% SET PATHs


% addpath ~/CancerGenomeAnalysis/trunk/matlab
% addpath ~/CancerGenomeAnalysis/trunk/matlab/snp
% addpath ~/CancerGenomeAnalysis/trunk/matlab/gp_modules



%% HANDLE ARGUMENTS

method_st = 'preprocess';

if isempty(varargin)
    disp(['Usage: gp_gistic_preprocessing -b base_dir -si sample_info_file ' ...
        '-al array_list_file -i input_snp_filename(s) -o output_snp_filename'...
        '-v verbose_level']);
    return
end
% a = handle_args({'b','si','al','ns','nc','ncn','up','bc', ...
%     'bems','bebpv','beapv','ss','svr','hn','ht','sh','cm',...
%     'i','o','a','nl','ec','h','sz','sm','mfn'},varargin);

a = handle_args({'o','up','bc','ht','ns','ncn','i','al','si','bems','beapv','v'},varargin);
% 
% if ~isempty(a.ss)
%     snp_skip=str2num(a.ss);
% else
     snp_skip=1;
% end
% 
% if ~isempty(a.sh)
%     show_hist = str2num(a.sh);
% else
     show_hist = 0;
% end

if ~isempty(a.v) && isnumeric(a.v)
    set_verbose_level(a.v);
end

array_list_file=a.al;
if isempty(array_list_file)
    error('must supply an array list file');
end


sample_info_file=a.si;
if isempty(sample_info_file)
    error('must supply a sample info file');
end

% 
% save_raw=a.svr;
% if isempty(save_raw)
     save_raw=0;
% else
%     save_raw=str2num(save_raw);
% end
% 
% norm_collapse_type=a.nc;
% if isempty(norm_collapse_type)
    norm_collapse_type=0; % mean
% end
norm_collapse_types={0,'mean';1,'median';2,'tangent'};
norm_collapse_method=enum_param(norm_collapse_type,norm_collapse_types);

norm_selection=a.ns;
if isempty(norm_selection)
    norm_selection=1; % closest_n
end
norm_selections={0,'all';1,'closest_n'};
norm_selection_method=enum_param(norm_selection,norm_selections);

use_paired=a.up;
if ~isempty(use_paired)
    use_paired=str2num(use_paired);
else
    use_paired=0; % do not used paired
end

if ~isempty(a.ncn)
    n_closest_n=str2num(a.ncn);
    if isempty(n_closest_n) || n_closest_n~=round(n_closest_n)
        error('n closest normals is not an integer number');
    end
else
    n_closest_n=5;
end

if ~isempty(a.bc)
    perform_batch_correction=str2num(a.bc);
else
    perform_batch_correction=0;
end

if ~isempty(a.bems)
    batch_effect_min_sz = str2num(a.bems);
else
    batch_effect_min_sz = 5;
end
% 
% if ~isempty(a.bebpv)
%     batch_effect_bonf_pv_thresh = str2num(a.bebpv);
% else
%     batch_effect_bonf_pv_thresh = .05;
% end

if ~isempty(a.beapv)
    batch_effect_abs_pv = str2num(a.beapv);
else
    batch_effect_abs_pv = .001;
end

base_dir=a.o;
base_dir=add_slash_if_needed(base_dir);
if ~exist(base_dir)
    mkdir(base_dir)
end

output_dir = add_slash_if_needed(base_dir);


% 
% if ~isempty(a.hn)
%     hist_qc_normals = str2num(a.hn);
% else
     hist_qc_normals = 0;
% end


if ~isempty(a.ht)
    hist_qc_tumors = str2num(a.ht);
else
    hist_qc_tumors = 1;
end


% 
% if ~isempty(a.cm)
%     conserve_memory = str2num(a.cm);
% else
     conserve_memory = 0;
% end


if isempty(a.i) || isempty(a.o)
    error('must have input and output files');    %Check that input and output files are specified
elseif ~isempty(regexp(a.i,'\.zip$', 'once' ))
    infile = unzip(a.i,base_dir);
else
    infile=regexp(a.i,'[^;]+','match');                                   %Assign filenames to variables INFILE and OUTFILE from information structure a
end
outdir=a.o;
[ofdir1,ofdir2] = fileattrib(outdir);
if ~ofdir1 || ~ofdir2.UserWrite
    error('Cannot write output snp file to specified directory')
end
%
% if ~isempty(a.a)
%     has_allele_data=str2num(a.a);                 %Assign whether allele data is provided
% else
     has_allele_data=0;
% end
% 
% if ~isempty(a.nl)                               %Assign number of lines to NLINES from A
%     nlines=str2num(a.nl);
%     if isempty(nlines)
%         error('number of lines is not a number');
%     end
% else
     nlines=-1;
% end
 
% if ~isempty(a.ec)                               %Assign number of extra columns to EXTRA_COLUMNS
%     extra_columns=str2num(a.ec);
%     if isempty(extra_columns) || extra_columns<0
%         error('extra columns should be a non-negative number');
%     end
% else
     extra_columns=0;
% end
 
% if ~isempty(a.h)
%     header_rows=str2num(a.h);                     %Assign number of header rows
%     if isempty(header_rows) || header_rows<0
%         error('header rows should be a non-negative number');
%     end
% else
     header_rows=0;
% end

% if ~isempty(a.sz)
%     read_chunks=str2num(a.sz);                    %Assign read chunks
%     if isempty(read_chunks) || read_chunks<0
%         error('read chunks should be a non-negative number');
%     end
% else
     read_chunks=10000;
% end


% if ~isempty(a.sm)
%     save_mat = a.sm;
% else
%     save_mat = 0;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALL SNP_TO_D

verbose(['Reading ' infile],10);

for k = 1:length(infile)

    M{k}=snp_to_D(infile{k},nlines,has_allele_data,extra_columns,header_rows,read_chunks);

end




%% CALL GISTIC_PREPROCESSING

[D,repnames,bad_tumor_names,arrays_for_core,celllines] = gistic_preprocessing(outdir,sample_info_file,array_list_file,M,...
    snp_skip,show_hist,0,norm_collapse_method,use_paired,...
    n_closest_n,perform_batch_correction,batch_effect_min_sz,[],...
    batch_effect_abs_pv,[],hist_qc_normals,hist_qc_tumors, conserve_memory,[],[],[],[],[],1);



s = whos('D');

if ~exist(output_dir,'dir')
    mkdir(output_dir);
end
% 
% dfile = [output_dir 'D.mat'];
% if save_mat
%     if s.bytes < 2000000000  %save without 'v7.3' tag is faster; '-v7.3' only necessary for files over 9 GB
%         save(dfile,'M');
%     else
%         save(dfile,'M','-v7.3');
%     end
% end

dupsfile = [output_dir 'duplicates'];
celllinesfile = [output_dir 'cell_lines'];
if ~isempty(repnames)
repnameschar = cellfun(@strvcat,repnames,'UniformOutput',0);
end

if ~isempty(repnames)
    fid = fopen(dupsfile,'w');
    for k = 1:length(repnames)
        fprintf(fid,'Appear to be duplicates: \n %s \n',repnameschar{k});
    end
    fclose(fid);
else
     fid = fopen(dupsfile,'w');
        fprintf(fid,'No Duplicates Found');
    fclose(fid);
end

% 
% if ~isempty(celllines) 
%     fid = fopen(celllinesfile,'w');
%     fprintf(fid,'Cell Lines: \n');
%     for k = celllines
%         fprintf(fid,'%s \n',D.sdesc{k});
%     end
%     fclose(fid);
% else
%       fid = fopen(celllinesfile,'w');
% 
%         fprintf('No Cell Lines Found');
% 
%     fclose(fid);
% end

badtumorsfile = [output_dir 'badtumors'];
fid = fopen(badtumorsfile,'w');
fprintf(fid,'Bad Tumors: \n');
fprintf(fid,'%s \n',bad_tumor_names{:});
fclose(fid);

S = cell2struct(arrays_for_core,'array',2);
arraylistout = [output_dir 'core_array_list_' datestr(now,'yymmdd') '.txt'];
infofilewrite(arraylistout,S);


%% WRITE SNP FILE(S)


cd(base_dir)
outfile = [output_dir 'preproc.snp'];
%write_as_dchip(outfile,D,1,[],2);
datwrite = mat2cell(D.dat,ones(1,size(D.dat,1)),ones(1,size(D.dat,2)));
poswrite = [D.marker cellfun(@num2str,mat2cell(D.chrn,ones(1,length(D.chrn)),1),'UniformOutput',0) ...
    cellfun(@num2str,mat2cell(D.pos,ones(1,length(D.pos)),1),'UniformOutput',0)];
if size(D.sdesc,1) > size(D.sdesc,2)
    D.sdesc = D.sdesc';
end


eol = '\r\n';  %for compatibility with pc operating systems
fid = fopen(outfile,'w');
verbose(['Writing output to:' outfile],20);
fprintf(fid,['%s\t%s\t%s' repmat('\t%s',1,length(D.sdesc)) eol],'SNP','Chromosome','PhysicalPosition',D.sdesc{:});

allwrite = [poswrite datwrite];
allwrite = reshape(allwrite',1,numel(allwrite));
formatstring = ['%s\t%s\t%s' repmat('\t%2.2f',1,length(D.sdesc)) eol];
fprintf(fid,formatstring,allwrite{:});
fclose(fid);



param_struct = struct('base_dir',base_dir,...
    'output_snp_file',outfile,'output_array_list_file',arraylistout,...
    'duplicates_file',dupsfile,'bad_tumors_file',badtumorsfile,...
    'sample_info_filename',sample_info_file, ...
    'array_list_filename',array_list_file, ...   
    'n_closest_n',num2str(n_closest_n), ...
    'use_paired_normal',num2str(use_paired),...
    'correct_batch_effect', ...
    num2str(perform_batch_correction),...
    'batch_effect_min_sz',num2str(batch_effect_min_sz),...
    'batch_effect_absolute_pv',num2str(batch_effect_abs_pv),...
    'hist_qc_tumors',num2str(hist_qc_tumors));
% 
%     'norm_collapse_method',norm_collapse_method,...
%      'batch_effect_bonf_pv',num2str(batch_effect_bonf_pv_thresh),...
%    'snp_skip',num2str(snp_skip),   
% 'save_raw', ...
%     num2str(save_raw),'hist_qc_normals',num2str(hist_qc_normals),...
% 'show_hist',num2str(show_hist));

param_file = [output_dir 'parameters.txt'];
verbose(['Writing parameter file to:' param_file],10);

gp_write_params(method_st,param_struct,[],param_file);








