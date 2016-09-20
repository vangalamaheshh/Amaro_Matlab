function regs=LOH_gistic(CNsort, extension, cnv_file, refgene, remove_X, r , focal_analysis, cap_vals, gen_regs_params)

if isempty(r)
    r=1 %leave 'r' out, default = 1
end

if ~isempty(cap_vals)
CNsort.dat(CNsort.dat<cap_vals(1))=cap_vals(1);
CNsort.dat(CNsort.dat>cap_vals(2))=cap_vals(2);
end

if remove_X
  CNsort=reorder_D_rows(CNsort,find(CNsort.chrn<=22));
end
CNsort=rmfield_if_exists(CNsort,{'cbs','cbs_rl'});

disp('Removing CNVs')
if exist('cnv_file','var') && ~isempty(cnv_file)
    CNsort=remove_cnv(CNsort,cnv_file);
end

base_dir='./';

if exist('refgene','var') && ~isempty(refgene)
    refgene_file=refgene;
else
    refgene_file  ='/xchip/tcga/Annotation_Files/UCSC_hg18/hg18_20070817.mat';
end

qv_thresh=0.25;
ext='.cneutral_loh';
t_amp = 0.1699;
t_del = 0.1926;
save_seg_data = 1;
res=0.0001; % (may change in the future default 0.001)
ts=[t_amp t_del];
method_st='.allele';

load(refgene_file);% read rg and cyto
rg=add_chrn(rg);
rg=mark_futreal_genes(rg);

if focal_analysis
  [mx,mn]=get_high_cutoffs(CNsort,cyto);
  fid=fopen([base_dir 'focal_sample_specific_cutoffs' ext '.txt'],'w');
  fprintf(fid,'Sample\tAmp cutoff\tDel cutoff\n');
  for i=1:size(CNsort.dat,2)
    fprintf(fid,'%s\t%f\t%f\n',CNsort.sdesc{i},mx(i),mn(i));
  end
  fclose(fid);
  CNsort.dat(CNsort.dat>0 & CNsort.dat<repmat(mx+ts(1),size(CNsort.dat,1),1))=0;
  CNsort.dat(CNsort.dat<0 & CNsort.dat>repmat(mn-ts(2),size(CNsort.dat,1),1))=0;
end

%% GISTIC CORE  modified
score_type=struct('method','nxa','amp_thresh',t_amp,'del_thresh',-t_del,'res',res);
%param_struct = run_gistic(base_dir,CNsort,refgenefile,array_list_file,t_amp,t_del,ext,qv_thresh,res,focal_analysis,gen_regs_params,genepattern);
%param_struct = run_gistic(base_dir,CNsort,refgene_file,[],t_amp,t_del,ext,qv_thresh,res,0,[],0);

[q,p,d,ads]=snp_score_permutations(CNsort,score_type,-1);
k=1;
if(isempty(find(q{k}<=qv_thresh, 1)))
    score_thresh(k)=max(ads{k})+0.01;%#ok
else
    score_thresh(k)=min(ads{k}(find(q{k}<=qv_thresh)));
end
score_thresh(2)=1;
CNsort=add_cyto(CNsort,cyto);

if isempty(gen_regs_params)
    regs=generate_regs_by_peel_off(CNsort,ads,d,q,score_type,score_thresh, struct('method','leave-k-out','k',r));
else
    arm_rl=runlength(CNsort.armn,CNsort.chrn);
    arm_rl=[arm_rl CNsort.chrn(arm_rl(:,1))];
    if isfield(gen_regs_params,'wide_type')
        regs=generate_regs_by_peel_off(CNsort,ads,d,q,score_type,score_thresh,...
            gen_regs_params.wide_type,[],...
            struct('method',gen_regs_params.peel_off_method,'regions_rl', ...
            arm_rl));
    else
        regs=generate_regs_by_peel_off(CNsort,ads,d,q,score_type,score_thresh,...
            struct('method','leave-k-out','k',gen_regs_params.k),[],...
            struct('method',gen_regs_params.peel_off_method,'regions_rl', ...
            arm_rl));
    end
%     regs=generate_regs_by_peel_off(CNsort,ads,d,q,score_type,score_thresh,...
%         struct('method','leave-k-out','k',gen_regs_params.k),[],...
%         struct('method',gen_regs_params.peel_off_method,'regions_rl', ...
%         arm_rl));
end
broad_type=struct('method','scorecutoff','ads',{ads},'p_arm',0.5,...
    'score_thresh',score_thresh,'score_thresh_focal',score_thresh);
%regs=find_broad_regs(CNsort,cyto,regs,broad_type);

%% sort regs by position
for k=1:2
  if ~isempty(regs{k})
      pos=cat(1,regs{k}.peak);
  [spos,sposi]=sort(pos);
  regs{k}=regs{k}(sposi);
  end
end
pvs=q;

%% output

core_out_file = [base_dir 'gistic_core' ext '.mat'];
chrn=CNsort.chrn; 
pos=CNsort.pos; 
save(core_out_file,'regs','ts','pvs','p','q','ads','chrn','pos');

%plot_snp_score([],[extension ext],CNsort,q,ads,0.25,1,1,[],cyto,[],[],[],[],[],0,regs,0,[],1); % 0
 make_all_lesions_file(['all_lesions' method_st ext extension '.txt'],CNsort,regs,pvs,cyto,[],ts,ts,1,0.5);
 all_lesions_file=['all_lesions' method_st ext extension '.txt'];
%gistic_plots(base_dir,plotsfname,CNsort,q,ads,regs,cyto,all_lesions_file,[],[],[],[],[],qv_thresh,[],[],[],[],genepattern); % qv_scale);
 gistic_plots(base_dir,[extension ext],CNsort,q,ads,regs,cyto,all_lesions_file,[],[],[],[],[],qv_thresh,[],[],[],[],1); % qv_scale);

seg_out_file = [base_dir 'segmented_data' ext '.mat'];
D=CNsort;
save(seg_out_file, 'D')
