function [param_struct,CL] = run_gistic(base_dir,CL,refgenefile,array_list_file,...
    t_amp,t_del,ext,qv_thresh,res,focal_analysis,gen_regs_params,genepattern)
% RUN_GISTIC -- Perform gistic CORE, gistic OUTFILES, and gistic PLOTS;
% (gp_gistic without the snp file conversion)
%
%[param_struct] = gistic(base_dir,C,segfile,refgenefile,...
%    cnv_file,array_list_file,t_amp,t_del,jss,ext,segheaderlines)
%
%
%base_dir: required
%CL (required = segmented data structure)
%refgenefile (required; see: '~gadgetz/projects/snp/data/Refgene/hg16/ucsc_20070112/hg16_20070112.mat')
%array_list_file (default: empty, use all)
%t_amp (default: .1)
%t_del (default: .1)
%ext (file extension; default = '')
%qv_thesh (default = 0.25);
%res (default = 0.001);
%genepattern (0) : outputs should be for genepattern module
%
% ---
% $Id$
% $Date: 2010-08-17 13:00:52 -0400 (Tue, 17 Aug 2010) $
% $LastChangedBy: stransky $
% $Rev$

%% Check inputs

varlist1 = {'base_dir','CL','refgenefile','array_list_file','t_amp','t_del','ext','qv_thresh','res','focal_analysis','gen_regs_params','genepattern'};

defaults = {'ERR','ERR','ERR','[]','0.1','0.1','''''','.25','0.001','0','[]','0'};

required = [1,1,1,0,0,0,0,0,0,0,0,0];

for idx = 1:length(varlist1)
    if required(idx) && (~exist(varlist1{idx},'var') || eval(['isempty(' varlist1{idx} ')']))
        error('Required input %s undefined.',varlist1{idx})
    elseif ~exist(varlist1{idx},'var') || eval(['isempty(' varlist1{idx} ')'])
        eval([varlist1{idx} '=' defaults{idx} ';'])
    end
end

base_dir = add_slash_if_needed(base_dir);

ts = [t_amp t_del];
all_lesions_file = [base_dir 'all_lesions_file' ext '.txt'];
scores_file = [base_dir 'scores'  ext '.gistic.txt'];


if ~iscell(CL.sdesc)
    CL.sdesc = cellstr(CL.sdesc);
end

load(refgenefile);  %loads cyto and rg

if focal_analysis
  [mx,mn]=get_high_cutoffs(CL,cyto);
  
  fid=fopen([base_dir 'focal_sample_specific_cutoffs' ext '.txt'],'w');
  fprintf(fid,'Sample\tAmp cutoff\tDel cutoff\n');
  for i=1:size(CL.dat,2)
    fprintf(fid,'%s\t%f\t%f\n',CL.sdesc{i},mx(i),mn(i));
  end
  fclose(fid);
  
  CL.dat(CL.dat>0 & CL.dat<repmat(mx+ts(1),size(CL.dat,1),1))=0;
  CL.dat(CL.dat<0 & CL.dat>repmat(mn-ts(2),size(CL.dat,1),1))=0;
end

%% match to array list file if exists

if exist('array_list_file','var') && ~isempty(array_list_file)
  AL = read_array_list_file(array_list_file);
  use_arrays = {AL.array};
  
  if isfield(AL,'inc_core')
    clear AL_core
    for i = 1:length(AL)
      AL_core(i)=str2num(AL(i).inc_core); %#ok<AGROW,ST2NM>
    end
    use_in_core=use_arrays(find(AL_core>0)); %#ok<FNDSB>
    use_in_outfiles=use_arrays(find(AL_core>-1)); %#ok<FNDSB>
  else
    use_in_core=use_arrays;
    use_in_outfiles=use_arrays;
  end
  [dum1,kpidx] = intersect(CL.sdesc,use_in_outfiles);
  CL = reorder_D_cols(CL,kpidx);
  if length(dum1) ~= length(use_in_outfiles)
    warning('Did not find a match to all arrays in array info file') %#ok<WNTAG>
  end
  [dum11,kpidx1] = intersect(CL.sdesc,use_in_core);
  CL1 = reorder_D_cols(CL,kpidx1);
  verbose(['Mathced ' num2str(length(kpidx1)) ' (out of ' num2str(size(CL.dat,2)) ') samples to array list file'],10);
else
  CL1 = CL;
  verbose('No array list file ... using all samples',10);
end

%% save data going into GISTIC core
if ~genepattern
    write_seg_file([base_dir 'input_to_gistic' ext '.seg.txt'],CL1,0,1);
    % save_D([base_dir 'input_to_gistic' ext '.mat'],CL1,'-v7.3');  %%% Removed due to problems with Broad and Stanford .chr field
end

%% GISTIC CORE  (1)

score_type=struct('method','nxa','amp_thresh',ts(1),'del_thresh',-ts(2),'res',res);

[q,p,d,ads]=snp_score_permutations(CL1,score_type,-1);

for k=1:2
  if(isempty(find(q{k}<=qv_thresh, 1)))
    score_thresh(k)=max(ads{k})+0.01;%#ok
  else
    score_thresh(k)=min(ads{k}(find(q{k}<=qv_thresh)));%#ok
  end
end
CL1=add_cyto(CL1,cyto);
if isempty(gen_regs_params)
  regs=generate_regs_by_peel_off(CL1,ads,d,q,score_type,score_thresh,...
                                 struct('method','leave-k-out','k',1));
  
  %% to use robust:
  %% regs=generate_regs_by_peel_off(CL1,ads,d,q,score_type,score_thresh,...
  %                               struct('method','robust','conf_level',[0.95
  %                               0.95],'nperm',10,'ampdel',[1 1]));
  %ampdel should be [1 0] for LOH.
  
else
  arm_rl=runlength(CL1.armn,CL1.chrn);
  arm_rl=[arm_rl CL1.chrn(arm_rl(:,1))];
  if isfield(gen_regs_params,'wide_type')
    regs=generate_regs_by_peel_off(CL1,ads,d,q,score_type,score_thresh,...
                                   gen_regs_params.wide_type,[],...
                                   struct('method',gen_regs_params.peel_off_method,'regions_rl', ...
                                          arm_rl));
  else
    regs=generate_regs_by_peel_off(CL1,ads,d,q,score_type,score_thresh,...
                                   struct('method','leave-k-out','k',gen_regs_params.k),[],...
                                   struct('method',gen_regs_params.peel_off_method,'regions_rl', ...
                                          arm_rl));
  end
end  
  % ,struct('method','regions','regions_rl',runlength(CL1.armn));


broad_type=struct('method','scorecutoff','ads',{ads},'p_arm',0.5,...
    'score_thresh',score_thresh,'score_thresh_focal',score_thresh);
regs=find_broad_regs(CL1,cyto,regs,broad_type);

%% sort regs by position
if isempty(regs{2})
	disp('Empty Deletion data.  Copying Amp data for compatibility.')
	regs{2} = regs{1};
	ads{2} = ads{1};
	broad_type.ads{2} = broad_type.ads{1};
	p{2} = p{1};
	q{2} = q{1};
end
for k=1:2
  pos=cat(1,regs{k}.peak);
  [spos,sposi]=sort(pos);
  regs{k}=regs{k}(sposi);
end


pvs=q;

core_out_file = [base_dir 'gistic_core' ext '.mat'];
chrn=CL1.chrn; %#ok<NASGU>
pos=CL1.pos; %#ok<NASGU>
if ~genepattern
    save(core_out_file,'regs','ts','pvs','p','q','ads','chrn','pos','d');
end
% should not save ext, it should be a parameter to GISTIC_outfiles --> Done
% added chrn and pos to allow to map to genome

clear CL1;

%% GISTIC OUTFILES (2)


%%%%%%%%%%%%%%%%%%%

%don't supply p=paramsfile, give ui and base_dir instead
% args_outfiles = struct('ui',CL ,'GISTIC','gistic_core.mat' ,'p',[],...   
%     'basedir',base_dir,'lout',all_lesions_file,'sout',scores_file ,'add_vals',1,...
%     'name', '','refgene_file',refgenefile,...
%     'annot_file',annotations_file,'mark',annotation_marker_symbol );
% this should get arguments, not a structure
% GISTIC_outfiles(args_outfiles);
if ~genepattern
    % write focal_broad.txt
    f=fopen([base_dir 'focal_broad' ext '.txt'],'w');
    ampdel={'Amp','Del'};
    for k=1:2
        for i=1:length(regs{k})
            [st,chr,bp]=genomic_location(CL,{[regs{k}(i).st regs{k}(i).en]},cyto,1,0); %#ok<NASGU>
            fprintf(f,'%s\n',[ ampdel{k} num2str(i) ')' st ' focal:' num2str(regs{k}(i).focal) ' broad:' num2str(regs{k}(i).broad)]);
        end
    end
    fclose(f);
end

make_all_lesions_file(all_lesions_file,CL,regs,pvs,cyto,[],ts,ts,1,broad_type);

verbose('Adding gene annotations',10);
rg=add_chrn(rg); %#ok<NODEF>
for i=1:length(rg)
    rg(i).symbol=rg(i).symb;
end
% annot_file=a.annot_file;
% mark=a.mark;
%rg=add_annotation(rg,annot_file,mark);
calls=call_regs(CL,regs,{ts}); 

verbose('Writing output gene table files',10);
genetables(CL,[],[],rg,cyto,regs,calls,ts,ext,all_lesions_file,base_dir,genepattern);

verbose('Writing q value, p value and amp/del score file',10);

if length(CL.pos)==length(q{1})
    sout=scores_file ;
    write_score_file(sout,CL,p,q,ads,ts);
else
    verbose(['Length of struct does not match length of scores!' num2str(size(CL.pos,2)) ' ' num2str(length(q{1}))],10);
    %    Score=[q,p,ads];
end

% copy_number_statistics(CL, ext, base_dir)

if ~genepattern
    write_bed_file([base_dir 'regions_track' ext '.bed'],'GISTIC','0,0,255',regs,CL,cyto);
end

%% GISTIC PLOTS (3)



plotsfname = [ '.' num2str(size(CL.dat,2)) ext ]; 

gistic_plots(base_dir,plotsfname,CL,q,ads,regs,cyto,all_lesions_file,[],[],[],[],[],qv_thresh,[],[],[],[],genepattern); % qv_scale);

 %% Make param structure

param_struct = struct('base_directory', base_dir,'qv_threshold',num2str(qv_thresh),'amplifications_threshold',num2str(t_amp),'deletions_threshold',num2str(t_del),...
    'array_list_file',array_list_file,'all_lesions_file',all_lesions_file,'scores_file',scores_file,'refgenefile',refgenefile,...
    'GISTIC_CORE_FILE',core_out_file,'res',res);
