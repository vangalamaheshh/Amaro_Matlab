function gp_gistic_preprocessing(varargin)
% gp_snp_subselection -b base_dir -dd data_dir -si sample_info_file -al array_list_file -o output_dir_ext
% -of output_file_name -ns norm_selection -nc norm_collapse_type -ncn n_closest_n -up use_paired_norm 
% -bc correct_batch_effect -bems batch_effect_min_sz -bebpv batch_effect_bonf_pv -beapv batch_effect_absolute_pv 
% -p parameter_file -ss snp_skip -svr save_raw
%
% input: a directory of .D.mat files (for the various plates) 
%        
%        a master sample_info file
%        a list of tumor and normal samples to use
%       

addpath /xchip/gistic/Code/GISTIC1.0/matlab/
addpath /xchip/gistic/Code/GISTIC1.0/matlab/snp
addpath /xchip/gistic/Code/GISTIC1.0/matlab/gp_modules

method_st='preprocess';
a=handle_args({'b','dd','si','al','o','of','ns','nc','ncn','up','bc','bems','bebpv','beapv','p','ss','svr'},varargin);

if ~isempty(a.ss)
  snp_skip=str2num(a.ss);
else
  snp_skip=1;
end

array_list_file=a.al;

sample_info_file=a.si;
if isempty(sample_info_file)
  error('must supply a sample info file');
end

datadir=a.dd;
if isempty(datadir)
  error('must supply a data directory');
end

save_raw=a.svr;
if isempty(save_raw)
  save_raw=0;
else
  save_raw=str2num(save_raw);
end

norm_collapse_type=a.nc;
if isempty(norm_collapse_type)
  norm_collapse_type=1; % median
end
norm_collapse_types={0,'mean';1,'median';2,'tangent'};
norm_collapse_method=enum_param(norm_collapse_type,norm_collapse_types);

norm_selection=a.ns;
if isempty(norm_selection)
  norm_selection=1; % closest_n
end
norm_selections={0,'all';1,'closest_n'};
norm_select_method=enum_param(norm_selection,norm_selections);

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
  n_closest_n=0;
end

if ~isempty(a.bc)
  perform_batch_correction=str2num(a.bc);
else
  perform_batch_correction=0;
end


base_dir=a.b;
if isempty(base_dir)
  base_dir='/';
else
  if base_dir(end)~='/'
    base_dir=[base_dir '/'];
  end
end

output_dir_extension=a.o;
if isempty(output_dir_extension)
  output_dir_extension='output';
end

%----- read sample info file
SI=read_sample_info_file([ sample_info_file '.txt']);
% add good field as yes if it does not exist
if ~isfield(SI,'good')
  for i=1:length(SI)
    SI(i).good='EMPTY';
  end
end

%----- read array list file
if ~isempty(a.al)
  use_arrays=read_dlm_file(array_list_file);
  use_arrays=cat(1,use_arrays{:});
  use_arrays=cellstr(unique(strvcat(use_arrays),'rows')); % make sure it is unique

  %----- match names to sample info
  [Mt,m1,m2]=match_string_sets_hash(use_arrays,{SI.array});
  if length(m1)<length(use_arrays) || isempty(m1) % did not match all the arrays
    error('Did not find a match in the sample info file to all the arrays in the array list file');
  end
  use_arrays_idx=m2;
else
  use_arrays_idx=1:length(SI);
end

% find files to load
if isfield(SI,'filename')
  plates=cellstr(unique(strvcat({SI(use_arrays_idx).filename}),'rows'));
else
  error('Sample info file must contain a column "filename"');
end

for i=1:length(plates)
  disp(plates{i});
  mat_file=[datadir plates{i} '.D.mat'];
  if exist(mat_file,'file')
    fprintf(1,['Reading matfile: ' mat_file newline]);
    x=load(mat_file);
    var_names=fieldnames(x);
    M{i}{1}=getfield(x,varnames{1}); % treat first variable in struct as D
    clear x
  else
    error(['Cannot find ' mat_file]);
  end
  [M{i}{1},IDX{i}{1}]=gistic_filter_samples(M{i}{1},SI,use_arrays_idx);
end
M=unite_Ds(cat(2,M{:}),'cols');

idx=cat(2,IDX{:});
idx=cat(1,idx{:});
SIs=SI(idx);  

%---- remove AFFY control probesets
affx_idx=grep('^AFFX',M.marker,1);
M1=reorder_D_rows(M,setdiff(1:length(M.marker),affx_idx));
M1=rmfield_if_exists(M1,'orig');

% subselect markers
snp_subset=1:snp_skip:size(M1.dat,1); 
M1=reorder_D_rows(M1,snp_subset);

% create a directory according to normalization type

output_dir=[ base_dir method_st '_' output_dir_extension '/'];
if ~exist(output_dir,'dir')
  mkdir(output_dir);
else
  error('Output directory already exists');
end

%----- match samples to sample_info and add genome_info
% FIXME: have a loop on the chip types (here only 1)
M2=add_sample_info(M,SIs,'array');

% [tmp,x1]=sortrows([ M2.supdat(find_supid(M2,'BATCH'),:)' M2.origidx']);
%[X,x1,x2]=match_string_sets(M2.sdesc,M1.sdesc);
%M2=reorder_D_cols(M2,x1);
M2=rmfield_if_exists(M2,'origidx');
M=M2;
clear M2

if save_raw
  save([ output_dir 'Mraw.D.mat'],M,'-v7.3');
end

%-- sort according to genomic location
M=order_by_pos(M);

params=struct('batch_effect_correction',struct('method','one_vs_all','min_sz','bonf_pv_threshold','absolute_pv'),...
              'collapse',norm_collapse_method,'selection',norm_select_method,...
              'use_paired',use_paired,'n_closest_n',n_closest_n,'perform_batch_correction',perform_batch_correction,...
              'use_all_cores_for_ref',1,'median_scale_before_batch_correct',1,'median_scale_after_norm',0);


[CL{1},P{1},batches{1},snps{1},mn{1}]=preprocess_snp_data(M,perform_batch_correction,1,0,n_close_normals,0,1,1, ...
                                                  use_tangent_normalization);


,remove_cnps);
[CL{2},P{2},batches{2},snps{2},mn{2}]=preprocess_snp_data(M,1,1,0,use_n_close_normals,0,1,1,use_tangent_normalization,remove_cnps);
[CL{1},P{1},batches{1},snps{1},mn{1}]=preprocess_snp_data(M,1,1,0,5); % put 1 at the end for matched normal
%-- remove markers that don't have a position (AFFX_*)
CL{1}=reorder_D_rows(CL{1},find(~isnan(CL{1}.pos)));
C=CL{1};


%-----------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------





use_n_close_normals=5;
use_tangent_normalization=0;
remove_cnps=0;

%---------- preprocess M (batch correction, divide by normals)
[CL{1},P{1},batches{1},snps{1},mn{1}]=preprocess_snp_data(reorder_D_cols(M,strmatch('2',{M.sis.ploidy})),1,1,0,use_n_close_normals,0,1,1,use_tangent_normalization,remove_cnps);
%-- remove markers that don't have a position (AFFX_*)
CL{1}=reorder_D_rows(CL{1},find(~isnan(CL{1}.pos)));
C=CL{1};

%-- sort according to genomic location
x=C.chrn*1e11+C.pos;
[sx,si]=sort(x);
C=reorder_D_rows(C,si);

clear CL1 CL

%---------- quality control (search for peaks in the histograms) 
throw_out_n=histogram_qc(C,['hist_qc_n.out.' method_st '.ps'],1);
throw_out_n=setdiff(throw_out_n,find(C.supdat(find_supid(C,'FORCE'),:)));

normals = strmatch({'2'},{C.sis.ploidy});
bad_normals = setdiff(normals, throw_out_n);
good_samples = setdiff(1:size(C.dat,2), bad_normals);

sdesc=C.sdesc;
sis=C.sis;

save sdesc_sis.mat sdesc sis
throw_out_n_names=C.sdesc(bad_normals);
save(['throw_out_n_normals.' method_st '.mat'],'bad_normals','throw_out_n_names');

if normal_hist_qc
normals = strmatch('2',{M.sis.ploidy});
M=reorder_D_cols(M,setdiff(1:size(M.dat,2),normals(bad_normals)));
end

clear C


%%%%%%%%%%%%%%%%%%%%%
%Second Time


%---------- preprocess M (batch correction, divide by normals)
[CL{1},P{1},batches{1},snps{1},mn{1}]=preprocess_snp_data(reorder_D_cols(M,strmatch('2',{M.sis.ploidy})),1,1,0,use_n_close_normals,0,1,1,use_tangent_normalization,remove_cnps);
%-- remove markers that don't have a position (AFFX_*)
CL{1}=reorder_D_rows(CL{1},find(~isnan(CL{1}.pos)));
C=CL{1};

%-- sort according to genomic location
x=C.chrn*1e11+C.pos;
[sx,si]=sort(x);
C=reorder_D_rows(C,si);

clear CL1 CL

%---------- quality control (search for peaks in the histograms) 
throw_out_n2=histogram_qc(C,['hist_qc_n2.out.' method_st '.ps'],1);
throw_out_n2=setdiff(throw_out_n,find(C.supdat(find_supid(C,'FORCE'),:)));

normals2 = strmatch({'2'},{C.sis.ploidy});
bad_normals2 = setdiff(normals2, throw_out_n2);
good_samples2 = setdiff(1:size(C.dat,2), bad_normals2);

sdesc2=C.sdesc;
sis2=C.sis;

save sdesc2_sis2.mat sdesc2 sis2
throw_out_n_names2=C.sdesc(bad_normals2);
save(['throw_out_n_normals2.' method_st '.mat'],'bad_normals2','throw_out_n_names2');

if normal_hist_qc
normals2 = strmatch('2',{M.sis.ploidy});
M=reorder_D_cols(M,setdiff(1:size(M.dat,2),normals2(bad_normals2)));  % craig 
end

clear C

%%%%%%%%%%%%%%%%%
%Third Time


[CL{2},P{2},batches{2},snps{2},mn{2}]=preprocess_snp_data(M,1,1,0,use_n_close_normals,0,1,1,use_tangent_normalization,remove_cnps);
%-- remove markers that don't have a position (AFFX_*)
CL{2}=reorder_D_rows(CL{2},find(~isnan(CL{2}.pos)));
CL{2}=reorder_D_cols(CL{2},strmatch(select_type,{CL{2}.sis.type}));
C=unite_Ds(CL);

%-- sort according to genomic location
x=C.chrn*1e11+C.pos;
[sx,si]=sort(x);
C=reorder_D_rows(C,si);

clear CL1 CL

throw_out=histogram_qc(C,['hist_qc.out.' method_st '.ps'],1);
throw_out=setdiff(throw_out,find(C.supdat(find_supid(C,'FORCE'),:)));
throw_out_names=C.sdesc(throw_out);
save(['throw_out.' method_st '.mat'],'throw_out','throw_out_names');

if tumor_hist_qc
C2=reorder_D_cols(C,setdiff(1:size(C.dat,2),throw_out));
else
C2=C;
end

clear C

cd([basedir dirname]);
save(['C2.' method_st '.mat'],'C2', '-v7.3');
copyqualityscore=calc_copy_quality(C,'absmedian');
save copyqualityscore.mat copyqualityscore
keyboard

