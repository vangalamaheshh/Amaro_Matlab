function gistic_pipeline_v2(select_rule,snpdir,data_fname,basedir,...
                            datadir,sample_info_dir,rundir,dirname,GI_fname,sample_info_fname,skip_snp,...
                            filter_snp,match_field_in_sample_info,refgene_file)

save([ basedir rundir 'start']);

close all
addpath ~gadgetz/matlab/snp

% basedir=['/xchip/projects/snp/'];
% dirname=['genentech/' select_type '/'];

use_cnag=0;
method_st='dChip';

%------- read snp files
if ischar(data_fname)
  clear tmp;
  tmp{1}=data_fname;
  data_fname=tmp;
end

for i=1:length(data_fname)
  no_dir_fname=regexprep(data_fname{i}(1:(end-3)),'.*/','');
  disp(no_dir_fname);
  mat_file=[basedir datadir no_dir_fname 'mat'];
  if exist(mat_file,'file')
    fprintf(1,['Reading ' mat_file newline]);
    x=load(mat_file);
    M{i}{1}=x.x;
  else
    fprintf(1,['Reading ' snpdir data_fname{i} newline]);
    if strcmp(data_fname{i}((end-2):end),'snp')
      has_chr_pos=1;
    else
      has_chr_pos=0;
    end
    M{i}{1}=read_modelled_data_file([ snpdir data_fname{i}],-1,1,has_chr_pos);
    x=M{i}{1};
    save(mat_file,'x');
  end
end
keyboard
N=unite_Ds(cat(2,M{:}),'cols');
M=N;
clear N; 
%---- we have M

%---- remove AFFY control probesets
M1=reorder_D_rows(M,75:size(M.dat,1));
M1=rmfield(M1,'orig');
M1=rmfield(M1,'affy_calls');
if ~exist('skip_snp','var') || isempty(skip_snp) || skip_snp==0
  skip_snp=1;
end
snp_subset=1:skip_snp:size(M1.dat,1); 
% snp_subset=1:size(M1.dat,1); 
M1=reorder_D_rows(M1,snp_subset);

%----- read genome_info file
% GI_fname=[ '~gadgetz/projects/snp/lung/Mapping250K_STY_genomeinfo_hg17_sorted'];
if exist([GI_fname '.mat'])
  load([GI_fname '.mat']);
else
  GI=read_genome_info_file([ GI_fname '.txt'],'fast');
  save([GI_fname '.mat'],'GI');
end

%----- read sample info file
SI_fname=[ basedir sample_info_dir sample_info_fname ];
SI=read_sample_info_file([ SI_fname '.txt']);

%----- select samples to use (from sample info)
if ischar(select_rule)
  clear tmp;
  tmp{1}={'+','type',select_type};
  select_type=tmp;
end

%dirname='';
g=[];
for i=1:size(select_rule,1)
  if select_rule{i,1}=='+'
    g=union(g,strmatch(lower(select_rule{i,3}),lower(strvcat(SI.(select_rule{i,2}))),'exact'));
  else
    g=setdiff(g,strmatch(lower(select_rule{i,3}),lower(strvcat(SI.(select_rule{i,2}))),'exact'));
  end
%  dirname=[dirname '_' select_type{i}];
end
% include rundir in dirname
%dirname=[ rundir dirname(2:end) '/' ]; % remove leading _
dirname=[ rundir dirname ]; 

if ~exist([ basedir dirname],'dir')
  mkdir([ basedir dirname]);
end

%- choose also the ploidy=2 samples (for normalization)
p2=findstrings(strvcat(SI(:).ploidy),'2')';

%- which are good samples?
if ~isfield(SI(1),'good')
  for i=1:length(SI)
    SI(i).good='yes';
  end
end
good=findstrings(strvcat(SI(:).good),'yes')';
% select all those that are g or p2 and are good
idx=unique(intersect(union(g,p2),good));

SIs=SI(idx);  
save([basedir dirname 'SIs.mat'],'SIs');
%-------------

cd([basedir dirname]);

%----- match samples to sample_info and add genome_info
% have a loop on the chip types (here only 1)
M2=add_snp_info(M1,GI,[],'by_name');    
M2=add_sample_info(M2,SIs,match_field_in_sample_info); %'array' or 'name'
[tmp,x1]=sort(M2.origidx);
%[X,x1,x2]=match_string_sets(M2.sdesc,M1.sdesc);
M2=reorder_D_cols(M2,x1);
M=M2;
clear M2
clear GI
%----------- M has all the information

keyboard

%---------- preprocess M (batch correction, divide by normals)
[CL{1},P{1},batches{1},snps{1},mn{1}]=preprocess_snp_data(M,1,1);
%-- remove markers that don't have a position (AFFX_*)
CL{1}=reorder_D_rows(CL{1},find(~isnan(CL{1}.pos)));
C=CL{1};

%-- remove samples which are marked as controls
not_controls=find(C.supdat(find_supid(C,'CTRL'),:)==0);
C=reorder_D_cols(C,not_controls);

%-- sort according to genomic location
x=C.chrn*1e11+C.pos;
[sx,si]=sort(x);
C=reorder_D_rows(C,si);
save C.mat C

clear CL1 M CL


%---------- quality control (search for peaks in the histograms)
throw_out=histogram_qc(C,['hist_qc.out.' method_st '.ps'],1);
throw_out=setdiff(throw_out,find(C.supdat(find_supid(C,'FORCE'),:)));
C2=reorder_D_cols(C,setdiff(1:size(C.dat,2),throw_out));
save(['throw_out.' method_st '.mat'],'throw_out');

C2=C;

% C2=C; % ignore qc
%--------- remove replicated samples (the bad copy)
C2=reorder_D_cols(C2,find(C2.supdat(find_supid(C2,'REP'),:)==0)); 

cd([basedir dirname]);
save(['C2.' method_st '.mat'],'C2');


%-----------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% smoothing: GLAD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use_samples_in=setdiff(1:length(C2.sdesc),resi);

glad_dir=[basedir dirname 'glad_' method_st];
if ~exist(glad_dir,'dir')
  mkdir(glad_dir);
end
cd(glad_dir);

run_glad_batch(C2,0,1,'normal',1);

% source run_all

% ----- merge all glad output files
f=fopen('output','w');
for i=1:size(C2.dat,2)
  seg=read_dlm_file(['Sample' sprintf('%03d',i) '.seg.dat']);
  for j=2:length(seg)
    fprintf(f,'%d %s %s %s %s %s\n',i,seg{j}{2},seg{j}{3},seg{j}{4},seg{j}{5},seg{j}{6});
  end
end
fclose(f);

if exist('use_samples_in','var')
  if ischar(use_samples_in)
    X=load(use_samples_in);
    [Mtmp,m1,idx]=match_string_sets(X.CL21.sdesc,C2.sdesc);
    if length(idx)~=size(X.CL21.dat,2)
      disp('no full match');
    end
  else
    idx=use_samples_in;
  end
else
  idx=[];
end

%cbs_dir=[basedir dirname 'cbs_' method_st];
%cd(cbs_dir);
[CL21,regs,pvs]=after_segmentation(C2,1,0,0,idx,filter_snp,[],[],[],[],['./output'],1);
write_as_dchip('segmented_data.txt',CL21,1);

%[CL21,regs,pvs]=after_segmentation(C2,1,0,0,idx,4,[],[],'_no_a',[],['./output'],1);
%make_all_lesions_file(['all_lesions.' method_st '.no_a.txt'],CL21,regs,pvs);



save CL21.mat CL21
% hdf5write('CL21.h5','/',CL21);

save regs.mat regs pvs
make_all_lesions_file(['all_lesions.' method_st '.txt'],CL21,regs,pvs);


if exist('refgene_file','var')
  load(refgene_file);
  % ~/projects/snp/data/Refgene/hg16/hg16.symb.mat
  % ~/projects/snp/data/Refgene/hg17/ucsc_20050723/hg17_20050723.mat
  rg=add_chrn(rg);
  
  futreal=read_dlm_file('~gadgetz/projects/snp/data/Futreal/Table_1_full_2006-02-16.txt');
  for i=2:length(futreal)
    fsymb{i-1}=futreal{i}{1};
    flocusid(i-1)=str2num(futreal{i}{3});
  end
  
  if ~isfield(rg,'symbol') && isfield(rg,'symb')
    for i=1:length(rg)
      rg(i).symbol=rg(i).symb;
    end
  end
    
%  [Mt,m1,m2]=match_string_sets_hash(fsymb,{rg.symb});
%  missing1=setdiff(1:length(fsymb),unique(m1));
  [Nt,n1,n2]=match_string_sets_hash(cellstr(num2str(flocusid')),cellstr(num2str([rg.locus_id]')));
  missing1a=setdiff(1:length(fsymb),unique(n1));
  un2=unique(n2);
  for i=1:length(un2)
    rg(un2(i)).symbol=[ rg(un2(i)).symbol '**'];
  end

  t=0.3;
  calls=call_regs(CL21,regs,t);
  tab=snp_generate_region_data(CL21,[],[],rg,regs,calls,t,pvs);

end

return



%--------------------------------
% LOH

loh_dir='/kinome/snp02/Barbara/dCHIPdata/500K/Sty/LOH_files/';
data_fname={'STY_ATLAS-ext_LOHdata.xls',...
            'STY_BIELD-ext_LOHdata.xls',...
            'STY_BLEAK-ext_LOHdata.xls',...
            'STY_BOTEL-ext_LOHdata.xls',...
            'STY_MYNAH_ext_LOHdata.xls',...
            'STY_OKAYS-ext_LOHdata.xls',...
            'STY_OSTIA-ext_LOHdata.xls',...
            'STY_REMEX-ext_LOHdata.xls',...
            'STY_STAND-ext_LOHdata.xls'};

data_mat_file='dChip_L_TSP_060608';
  
for i=1:length(data_fname)
  disp(i);
  if exist([basedir dirname data_fname{i}(1:(end-3)) 'loh.mat'],'file')
    fprintf(1,['Reading ' basedir dirname data_fname{i}(1:(end-3)) 'mat' newline]);
    x=load([basedir dirname data_fname{i}(1:(end-3)) '.loh.mat']);
    L{i}{1}=x.x;
  else
    fprintf(1,['Reading ' loh_dir data_fname{i} newline]);
    L{i}{1}=read_copy_number_file([ loh_dir data_fname{i}]);
    x=L{i}{1};
    save([basedir dirname data_fname{i}(1:(end-3)) '.loh.mat'],'x');
  end
end
% for i=1:length(data_fnames)
%      x=M{i}{1};
%      save([basedir dirname data_fname{i}(1:(end-3)) 'mat'],'x');
%  end
%  hdf5write( [pwd data_mat_file '.h5'],'/',M);
%  save([ basedir data_mat_file '.mat'],'M'); 
N=unite_Ds(cat(2,L{:}),'cols');
L=N;
clear N; 

clear all
close all
addpath ~/matlab/snp
basedir=['~/projects/snp/'];
dirname='lung/';
method_st='dChip';
glad_dir=[basedir dirname 'glad_' method_st];

cd(glad_dir);
load CL21
load regs

%% LOH
L=read_copy_number_file([ ' ' '100K glioma LOH_060419.txt']);
L.dat=L.dat(:,1:(end-1));
if same(L.marker(1:end-1),CL21.marker(1:end)) 
  L.marker=L.marker(1:(end-1));
  L.chr=L.chr(1:(end-1));
  L.pos=L.pos(1:(end-1));
  L.cM=L.cM(1:(end-1));
else
  disp('DO NOT MATCH');
end

[Mt,m1,m2]=match_string_sets(L.sdesc,CL21.sdesc);
L2=reorder_D_cols(L,m1);
L2=add_chrn(L2);

cd('LOH');

% discretize 
L2t=L2;
L2t.dat=double(L2t.dat>0.5);

score_type=struct('method','nxa','amp_thresh',0.3,'del_thresh',-0.3,'res',0.0001);
[q,p,d,ads]=snp_score_permutations(L2t,score_type,-1);
plot_snp_score(num2str(size(L2t.dat,2)),L2t,q,ads,0.25,1,1); % 0
qv_thresh=0.25;
k=1;
score_thresh(k)=min(ads{k}(find(q{k}<=qv_thresh)));
score_thresh(2)=1;
regs=generate_regs_by_peel_off(L2t,d,q,score_type,score_thresh,501);
pvs=q;

save L2.mat L2 L2t
save LOHregs.mat regs pvs

make_all_lesions_file(['all_lesions.' method_st '.LOH.txt'],L2,regs,pvs);

%-------------------------------------------------------------------------------------

[Mt,m1,m2]=match_string_sets(gcm,lst);
setdiff(1:length(lst),m2)

load('../../GCM/C2.dChip.mat');
C2=reorder_D_cols(C2,m1);
save(['C2.' method_st '.mat'],'C2');

for i=1:length(m2)
  unix_st=['ln -s ../../GCM/glad_dChip/' lst{m2(i)} '.dat .'];
  disp(unix_st);
  unix(unix_st);

  unix_st=['ln -s ' lst{m2(i)} '.dat Sample' sprintf('%03d',i) '.dat'];
  disp(unix_st);
  unix(unix_st);

  unix_st=['ln -s ../../GCM/glad_dChip/Sample' sprintf('%03d',m1(i)) '.seg.dat Sample' sprintf('%03d',i) '.seg.dat'];
  disp(unix_st);
  unix(unix_st);
end

unixst={'foreach fl (*.seg.dat)',...
        'tail +2 $fl >> output',...
        'end',...
        'awk ''BEGIN { x=23 } { if (($2==1) && (x==23)) { s+=1;} x=$2; print s,$2,$3,$4,$5,$6 }'' output > output.1',... 
        'mv output.1 output'};

disp('Execute this:');
disp(strvcat(unixst));

idx=[];
[CL21,regs,pvs]=after_segmentation(C2,1,0,0,idx,4,[],[],[],[],['./output'],1);
write_as_dchip('segmented_data.txt',CL21,1);

save CL21.mat CL21
save regs.mat regs pvs

make_all_lesions_file(['all_lesions.' method_st '.txt'],CL21,regs,pvs);
