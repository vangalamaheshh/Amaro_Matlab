function CL = struct_from_seg(base_dir,C,seg_dir,cnv_file,array_list_file,ext,segheaderlines,save_seg_data,remove_X,join_segment_size)
% struct_from_seg -- Perform gistic CORE, gistic OUTFILES, and gistic PLOTS;
% (gp_gistic without the snp file conversion)
%
%[param_struct] = gistic(base_dir,C,segfile,refgenefile,...
%    cnv_file,array_list_file,t_amp,t_del,jss,ext,segheaderlines)
%
%
%base_dir: required
%C (required = preprocessing output structure)
%seg_dir (required; directory with segmented data)
%cnv_file  (default: empty)  see:'/xchip/data01/gadgetz/CNV'
%array_list_file (default: empty, use all)
%ext (file extension; default = '')
%segheaderlines (default = 0);
%remove_X (default = 0)
%
% ---
% $Id$
% $Date: 2010-02-15 12:43:26 -0500 (Mon, 15 Feb 2010) $
% $LastChangedBy: rameen $
% $Rev$

%% Check inputs

varlist1 = {'base_dir','C','seg_dir','cnv_file','array_list_file','ext','segheaderlines','save_seg_data','remove_X','join_segment_size'};

defaults = {'ERR','ERR','ERR','[]','[]','''','0','1','0','1'};

required = [1,1,1,0,0,0,0,0,0,0];

for idx = 1:length(varlist1)
    if required(idx) && (~exist(varlist1{idx},'var') || eval(['isempty(' varlist1{idx} ')']))
        error('Required input %s undefined.',varlist1{idx})
    elseif ~exist(varlist1{idx},'var') || eval(['isempty(' varlist1{idx} ')'])
        eval([varlist1{idx} '=' defaults{idx} ';'])
    end
end

%% Convert C from HDF5
%if strcmp(class(C),'datastruct')
%  C_nohdf5 = convert_to_struct(C);
%else
%  C_nohdf5 = C;
%end
%
%clear C

%% Generate merged segmented data

%n = combine_segmented_data(C_nohdf5,seg_dir,'Sample');
n = combine_segmented_data(C,seg_dir,'Sample','.seg.dat',[seg_dir 'merged_Sample.seg.dat'],segheaderlines);

%% Generate CL

CL = make_D_from_seg([seg_dir 'merged_Sample.seg.dat'],[seg_dir 'GenCoords.dat']);
%CL = read_cbs_file([seg_dir 'merged_Sample.seg.dat'],C_nohdf5,segheaderlines,[],1,0,0);
%CL = read_cbs_file(segfile,C,segheaderlines,[],0,0,1);
%CL.dat = CL.cbs;

%%% Add .sis field

[M mi mj] = match_string_sets_hash(CL.sdesc,C.sdesc);

C = reorder_D_cols(C,mj);
CL = reorder_D_cols(CL,mi);

if isfield(C,'sis')
  CL.sis = C.sis;
end

clear C

%% Remove CNV
if exist('cnv_file','var') && ~isempty(cnv_file)
    CL=remove_cnv(CL,cnv_file);
end

%% Join small segments
if join_segment_size > 1
  CL.cbs=CL.dat;
  CL=smooth_cbs(CL,join_segment_size);
  CL.dat=CL.cbs;
  CL=rmfield_if_exists(CL,{'cbs','cbs_rl'});
end
  
%% remove X,Y chromosomes
if remove_X
  CL=reorder_D_rows(CL,find(CL.chrn<=22));
end
CL=rmfield_if_exists(CL,{'cbs','cbs_rl'});

%% subtract median
CL.medians=median(CL.dat(find(CL.chrn<=22),:),1);
CL.dat=CL.dat-repmat(CL.medians,size(CL.dat,1),1);

%% save segmented data
if save_seg_data
  seg_data_file = [base_dir 'segmented_data' ext '.mat'];
  save_D(seg_data_file,CL,'-v7.3');
end
