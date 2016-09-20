function gp_gistic(varargin)
%
% -b base_dir -sgfl segfile  [-ta t_amp] [-td t_del]
% -snpfile snpfile -refgen refgenefile [-alldat has_allele_data] 
%[-excols snpfile_has_extra_columns] [-snphead snp_header_rows] 
% [-qvt qv_thresh] [-alf array_list_file] [-js join_segment_size(def=5)] 
% [-cnv cnv_file] [-ext extension]
%---
% $Id$
% $Date: 2008-04-11 12:57:46 -0400 (Fri, 11 Apr 2008) $
% $LastChangedBy: jdobson $
% $Rev$

% addpath ~/CancerGenomeAnalysis/trunk/matlab/
% addpath ~/CancerGenomeAnalysis/trunk/matlab/gp_modules
% addpath ~/CancerGenomeAnalysis/trunk/matlab/snp

if isempty(varargin)
    disp('Usage: gp_gistic -b base_dir -sgfl segmentation_file ')
    disp('[-ta amplifications_threshold(def=.1)] [-td deletions_threshold(def=.1)] -snpfile input_snp_file')
    disp('-refgene refgenefile [-alldat snpfile_has_allele_data(def=0)]')
    disp('[-excols snpfile_extra_columns(def=0)] [-snphead snp_header_rows(def=0)]]')
    disp('[-qvt qv_thresh(def=0.25)] [-alf array_list_file(def:empty)] [-js join_segment_size(def=5)]')
    disp('[-cnv cnv_file] [-ext extension]');
    return
end

method_st = 'GISTIC';

a = handle_args({'b','sgfl','ta','td','snpfile','refgene','alldat',...
    'excols','snphead','qvt','alf','js','cnv','ext'},varargin);

base_dir = a.b;
if isempty(base_dir)
    error('Must supply a base directory')
end

%%%SEGMENTATION FILE PARAMETERS

segfile = a.sgfl;
if isempty(segfile)
    error('Must supply a segmentation file.  Use tag -sgfl.');
end

segheaderlines = 0;

%%%


if isempty(a.qvt)
    qv_thresh = 0.25;
else
    qv_thresh = str2num(a.qvt);
end


%Amplification and deletion thresholds
if isempty(a.ta)
    t_amp = .1;
else
    t_amp = str2num(a.ta);
end

if isempty(a.td)
    t_del = .1;
else
    t_del = str2num(a.td);
end


%%%SNP FILE PARAMETERS%%%

if isempty(a.snpfile)
%    error('Must supply a snpfile.  Use tag -snpfile.');
    snpfile = [];
else
    snpfile = a.snpfile;
end

if isempty(a.alldat)
    has_allele_data = 0;
else
    has_allele_data = a.alldat;
end

if ~isempty(a.excols)                               %Assign number of extra columns to EXTRA_COLUMNS
    extra_columns=str2num(a.excols);
    if isempty(extra_columns) || extra_columns<0
        error('extra columns should be a non-negative number');
    end
else
    extra_columns=0;
end

if ~isempty(a.snphead)
    header_rows=str2num(a.snphead);                     %Assign number of header rows
    if isempty(header_rows) || header_rows<0
        error('header rows should be a non-negative number');
    end
else
    header_rows=0;
end

%%% ARRAY LIST PARAMETERS

if ~isempty(a.alf)
    array_list_file = a.alf;
end

%%% JOIN SEGMENT SIZE
if isempty(a.js)
    join_segment_size = 5;
else
    join_segment_size = str2num(a.js);
end

%%% JOIN SEGMENT SIZE
if isempty(a.cnv)
    cnv_file='';
else
    cnv_file=a.cnv;
end

%%% JOIN SEGMENT SIZE
if isempty(a.ext)
    ext='';
else
    ext=a.ext;
end

%%% REFGENE FILE

if isempty(a.refgene)
    error('Must supply a refgene file.  Ex: -refgene hg17_20070131.mat');
else
    refgenefile = a.refgene;
end

%% READ Craw

verbose([datestr(now) ': Passing ' snpfile ' to read_modelled_data_file.m'],10);
Craw=read_modelled_data_file(snpfile,-1,-1,1,0,0,0,0,10000);

%%% DEFINE FILENAME EXTENSION

ext=['.' datestr(now,'yymmdd') ext];


%%% DON'T WRITE SEGMENTED DATA
save_seg_data = 0;


param_struct = run_gistic_from_seg(base_dir,Craw,segfile,refgenefile,cnv_file,array_list_file,t_amp,t_del,join_segment_size,ext,segheaderlines,qv_thresh,save_seg_data)





% CL = read_cbs_file(segfile,Craw,segheaderlines,[],1,0,1);
% CL.dat=CL.cbs;

% %% Remove CNV
% if exist('cnv_file','var') && ~isempty(cnv_file)
%     CL=remove_cnv(CL,cnv_file);
% end


% %% Join small segments
% CL.cbs=CL.dat;
% CL=smooth_cbs(CL,join_segment_size);
% CL.dat=CL.cbs;
% CL=rmfield_if_exists(CL,{'cbs','cbs_rl'});

% %% remove X,Y chromosomes

% CL=reorder_D_rows(CL,find(CL.chrn<=22));
% CL=rmfield_if_exists(CL,{'cbs','cbs_rl'});

% %% subtract median
% CL.medians=median(CL.dat,1);
% CL.dat=CL.dat-repmat(CL.medians,size(CL.dat,1),1);

% %% match to array list file if exists

% if exist('array_list_file','var')
%     AL = read_array_list_file(array_list_file);
%     use_arrays = {AL.array};
%     [dum1,kpidx,dum2] = intersect(CL.sdesc,use_arrays);
%     CL = reorder_D_cols(CL,kpidx);
%     if length(dum1) ~= length(use_arrays)
%         warning('Did not find a match to all arrays in array info file')
%     end
% end


% %% GISTIC CORE  (1)


% score_type=struct('method','nxa','amp_thresh',ts(1),'del_thresh',-ts(2), ...
%                   'res_per_sample',0.1,);

% [q,p,d,ads]=snp_score_permutations(CL,score_type,-1);

% for k=1:2
%     score_thresh(k)=min(ads{k}(find(q{k}<=qv_thresh)));
% end
% regs=generate_regs_by_peel_off(CL,ads,d,q,score_type,score_thresh,...
%     struct('method','leave-k-out','k',1));

% load(refgenefile);  %loads cyto and rg

% broad_type=struct('method','scorecutoff','ads',{ads},'p_arm',0.5,...
%     'score_thresh',score_thresh,'score_thresh_focal',score_thresh);
% regs=find_broad_regs(CL,cyto,regs,broad_type);

% pvs=q;


% core_out_file = [base_dir 'gistic_core' ext '.mat'];
% save(core_out_file,'regs','ts','pvs','p','q','ads','ext');
% % should not save ext, it should be a parameter to GISTIC_outfiles

% %% GISTIC OUTFILES (2)


% %%%%%%%%%%%%%%%%%%%
% annotation_marker_symbol='';
% annotations_file='';
% %don't supply p=paramsfile, give ui and base_dir instead
% % args_outfiles = struct('ui',CL ,'GISTIC','gistic_core.mat' ,'p',[],...   
% %     'basedir',base_dir,'lout',all_lesions_file,'sout',scores_file ,'add_vals',1,...
% %     'name', '','refgene_file',refgenefile,...
% %     'annot_file',annotations_file,'mark',annotation_marker_symbol );
% % this should get arguments, not a structure
% % GISTIC_outfiles(args_outfiles);

% % write focal_broad.txt
% f=fopen([base_dir 'focal_broad' ext '.txt'],'w');
% ampdel={'Amp','Del'};
% for k=1:2
%     for i=1:length(regs{k})
%         [st,chr,bp]=genomic_location(CL,{[regs{k}(i).st regs{k}(i).en]},cyto,1,0);
%         fprintf(f,'%s\n',[ ampdel{k} num2str(i) ')' st ' focal:' num2str(regs{k}(i).focal) ' broad:' num2str(regs{k}(i).broad)]);
%     end
% end
% fclose(f);

% make_all_lesions_file(all_lesions_file,CL,regs,pvs,cyto,[],ts,ts,1,broad_type);

% disp('Adding gene annotations');
% rg=add_chrn(rg);
% for i=1:length(rg)
%     rg(i).symbol=rg(i).symb;
% end
% % annot_file=a.annot_file;
% % mark=a.mark;
% % rg=add_annotation(rg,annot_file,mark);
% calls=call_regs(CL,regs,{ts});

% disp('Writing output gene table files');
% genetables(CL,[],[],rg,cyto,regs,calls,ts,ext,all_lesions_file,base_dir);

% disp('Writing q value, p value and amp/del score file');

% if length(CL.pos)==length(q{1})
%     sout=[scores_file ];
%     write_score_file(sout,CL,p,q,ads,ts);
% else
%     disp(['Length of struct does not match length of scores!' num2str(size(CL.pos,2)) ' ' num2str(length(q{1}))]);
%     %    Score=[q,p,ads];
% end
% copy_number_statistics(CL, ext, base_dir)



% %% GISTIC PLOTS (3)



% plotsfname = [ '.' num2str(size(CL.dat,2)) ext ]; 
% qv_scale=[];
% gistic_plots(base_dir,plotsfname,CL,q,ads,regs,cyto,all_lesions_file); % qv_scale);
% %% WRITE PARAMS

% %MODIFY THIS PART TO INCLUDE ALL METHODS


% param_struct = struct('base_directory', base_dir,'segmentation_file',segfile,'snp_file',snpfile,...
%     'qv_threshold',num2str(qv_thresh),'amplifications_threshold',num2str(t_amp),'deletions_threshold',num2str(t_del),...
%     'array_list_file',array_list_file,'all_lesions_file',all_lesions_file,'scores_file',scores_file,'refgenefile',refgenefile,...
%     'plots_file',plotsfname,'GISTIC_CORE_FILE',core_out_file);

param_file = [base_dir 'parameters_' method_st ext '.txt'];

disp(['Writing parameter file to:' param_file]);
local_copy=0;
gp_write_params(method_st,param_struct,[],param_file,local_copy);




