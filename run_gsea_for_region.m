function [G,m]=run_gsea_for_region(glad_dir,alg_dir,subset_dir,C,regs,pvs,rg,G,ad,regid,clsids,chr,lesion_name,nms,phens,no_gsea,force_run,t_vals)
% run GSEA for each region

addpath ~/matlab/gp_modules
addpath ~/matlab/snp
% G=add_D_sup(G,C.supacc,C.supdesc,C.supdat,'cols');
add_supid=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(clsids)<=1
    two_class=1;
    cls1id=clsids;
    if ~isempty(cls1id)
        cls1=C.supdat(cls1id,:);
    else
        cls1=zeros(1,size(C.dat,2));
    end
    cls2=[];
    cls0=cls1;
    cls0(~isnan(cls0))=1-cls0(~isnan(cls0));
else
    two_class=0;
    cls1id=clsids(1);
    cls2id=clsids(2);
    if ~isempty(cls2id)
        cls2=C.supdat(cls2id,:);
    else
        cls2=zeros(1,size(C.dat,2));
    end
    if ~isempty(cls1id)
        cls1=C.supdat(cls1id,:);
    else
        cls1=zeros(1,size(C.dat,2));
    end
    nanpos=find(isnan(cls1) | isnan(cls2));
    cls1(nanpos)=NaN;
    cls2(nanpos)=NaN;
    intersect_pos=find(cls1==1 & cls2==1);
    cls1(intersect_pos)=NaN;
    cls2(intersect_pos)=NaN;
    cls0=cls1+cls2;
    cls0(~isnan(cls0))=~cls0(~isnan(cls0));
end

if two_class
    [dum,idx]=max([ cls0; cls1],[],1);
    idx(find(isnan(cls0)))=NaN;
    [G,cur_supid]=add_D_sup(G,[ lesion_name ': 1-No/2-' nms{1} ],[lesion_name ': 1-No/2-' nms{1}],idx,'cols');
else
    [dum,idx]=max([ cls0; cls1; cls2],[],1);
    idx(find(isnan(cls0)))=NaN;
    [G,cur_supid]=add_D_sup(G,[ lesion_name ': 1-No/2-' nms{1} '/3-' nms{2}],[ lesion_name ': 1-No/2-' nms{1} '/3-' nms{2} ],idx,'cols');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove samples with NaN in any of the phenotypes
not_nan=find(~isnan(dum));
n_with_nan=length(find(isnan(dum)));
if n_with_nan>0
    disp(['removing ' num2str(n_with_nan) ' samples due to NaN in phenotype']);
end
[ss,si]=sort(idx(not_nan));
si=not_nan(si);
G=reorder_D_cols(G,si);
C=reorder_D_cols(C,si);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove the genes on the chromosome of the lesion
for ci=1:length(clsids)
  chrs{ci}=regexp(G.supdesc(clsids(ci),:),'chr([0-9XY]+)','tokens');
  chrs_num{ci}=chromosome2num(cat(1,chrs{ci}{:}));
end
chrs_num=unique(cat(1,chrs_num{:}));
chrs_num=chrs_num(~isnan(chrs_num));
if ~isempty(chrs_num)
  disp(['removing chromosomes ' num2str(chrs_num')]);
  not_chr=find(~ismember(G.chrn,[-1 chrs_num']));
  X=reorder_D_rows(G,not_chr);
else
  disp('Not removing any region');
  X=G;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cd([ glad_dir '/with_GE/' alg_dir '/' subset_dir ]);


tmp=sprintf('%d_',chrs_num);
chrst=tmp(1:(end-1));
dir_name=[ 'no' chrst];
if ~exist(dir_name,'dir')
    mkdir(dir_name);
end
cd(dir_name);

if ~exist(lesion_name,'dir')
    mkdir(lesion_name);
end
cd(lesion_name);

exp_fname=['gene_expression.no' chrst '.' lesion_name '.' num2str(size(X.dat,2)) '.gct'];
if ~exist(exp_fname,'file')
    %  write_mit_gct_file(exp_fname,X);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

if ~two_class
    options={['No_versus_' nms{1}],1,2,3;['No_versus_' nms{2}],1,3,2};
    for i=1:size(options,1)
        cls0=options{i,2};
        cls1=options{i,3};
        tst=options{i,4};

        n_cls0=length(find(X.supdat(cur_supid,:)==cls0));
        n_cls1=length(find(X.supdat(cur_supid,:)==cls1));
        n_tst=length(find(X.supdat(cur_supid,:)==tst));

        if (n_cls0 <= 2) || (n_cls1 <= 2)
            disp('One of the classes is too small');
            continue
        end

        
        v=X.supdat(cur_supid,:);
        [Dtrain,train_supid]=add_D_sup(X,'TRAIN','TRAIN',v==cls1,'cols');
        Dtrain=reorder_D_cols(Dtrain,find(ismember(v,[cls0 cls1])));

        Dtest=reorder_D_cols(X,find(ismember(X.supdat(cur_supid,:),[tst])));
        [Dtest,test_supid]=add_D_sup(Dtest,'TEST','TEST',ones(1,size(Dtest.dat,2)),'cols');

        if (0)

            [m.idx,m.q,m.p,m.s,m.pi0,m.F]=get_top_markers(Dtrain,train_supid,...
                struct('method','ttest'),...
                -1,...
                struct('method','bhfdr','thresh',0.1));
            %                                             struct('method','ttest_minvar','minvar',(0.3)^2),...
            %                                             struct('method','top','n',30));
            Xm=reorder_D_sup(reorder_D_rows(X,m.idx),'cols',[ find_supid(X,phenqs) cur_supid]);

            if length(m.idx)>5

                [Xm_ord,Xm_ord.gdend,Xm_ord.sdend]=two_way_clustering(Xm,...
                    struct('dist','correlation','preproc',...
                    'row_center_and_normalize','cluster','average'),...
                    struct('dist','euclidean','preproc',...
                    'row_center_and_normalize','cluster','ward'));



                display_D(Xm_ord,Xm_ord.gdend,Xm_ord.sdend,'dends');
                print_D([ lesion_name '_' options{i,1} '.FDR_0.1.clustering'],{{'pdf'}},1);
                disp([ lesion_name ' ' options{i,1} ' FDR 0.1: ' num2str(length(m.idx)) ' genes']);
            else
                disp([ lesion_name ' ' options{i,1} ' FDR 0.1: ' num2str(length(m.idx)) ' genes < 5']);
            end
        end
        if (0)
            trn_name=[ lesion_name '_' options{i,1} '.train'];
            tst_name=[ lesion_name '_' options{i,1} '.test'];

            if ~exist([ trn_name '.gct'],'file') || ~exist([ trn_name '.cls'],'file')
                write_mit_gct_file([ trn_name '.gct'],Dtrain);
                write_mit_cls_file([ trn_name '.cls'],Dtrain,train_supid,1,1);
            else
                disp(['Using existing ' trn_name '.gct']);
            end
            if ~exist([ tst_name '.gct'],'file') || ~exist([ tst_name '.cls'],'file')
                write_mit_gct_file([ tst_name '.gct'],Dtest);
                write_mit_cls_file([ tst_name '.cls'],Dtest,test_supid,1,1);
            else
                disp(['Using existing ' trn_name '.gct']);
            end

            if ~exist([ trn_name '.xvres.odf'],'file')
                %        gp_pnn_xvalidation([ trn_name '.gct'],[ trn_name '.cls'],...
                %                           [ trn_name '.pred.xv'],[ trn_name '.xvres.odf'],...
                %                           '2','0.3','10:10:300','2','1:0.5:4','1','0.5','1',[ trn_name '.xvres.mat']);
                gp_pnn_xvalidation([ trn_name '.gct'],[ trn_name '.cls'],...
                    [ trn_name '.pred.xv'],[ trn_name '.xvres.odf'],...
                    '2','0.3','10:10:300','2','1:0.5:4','1','0.5','0',[ trn_name '.xvres.mat']);
            else
                disp(['Using existing ' trn_name '.xvres.odf']);
            end

            if ~exist([ tst_name '.full.odf'],'file')
                gp_pnn([ trn_name '.gct'],[ trn_name '.cls'],...
                    [ tst_name '.gct'],[ tst_name '.cls'],[ tst_name '.full.odf'],...
                    [ tst_name '.pred.odf'],...
                    '10','2','0.3','2','1','1','0.5','euclidean','0',[ trn_name '.xvres.odf']);
            else
                disp(['Using existing ' tst_name '.full.odf']);
            end


            xv=read_mit_odf_file([ trn_name '.xvres.odf']);
            tst_res=read_mit_odf_file([ tst_name '.full.odf']);
            fet_nms=regexp(tst_res.data{1}{2},'(\w*)','tokens');
            fet_nms=cat(1,fet_nms{:});

            [Mt,m1,m2]=match_string_sets_hash(X.gacc,fet_nms);

            Xm=reorder_D_sup(reorder_D_rows(X,m1),'cols',[ find_supid(X,phens) cur_supid]);

            [Xm_ord,Xm_ord.gdend,Xm_ord.sdend]=two_way_clustering(Xm,...
                struct('dist','correlation','preproc',...
                'row_center_and_normalize','cluster','average'),...
                struct('dist','euclidean','preproc',...
                'row_center_and_normalize','cluster','ward'));


            display_D(Xm_ord,Xm_ord.gdend,Xm_ord.sdend,'dends');
            print_D([ lesion_name '_' options{i,1} '.xv.clustering'],{{'pdf'}},1);
            disp([ lesion_name ' ' options{i,1} ' xv: ' num2str(length(m1)) ' genes']);

            if (0)
              [confmat,nerr,cvres,D,params]=crossvalidate(reorder_D_cols(Xm,find(~ismember(Xm.supdat(end,:),tst))),...
                                                          size(Xm.supdat,1),...
                                                          struct('method','loocv','classifier',...
                                                                struct('method','matlab_classify',...
                                                                'dist_type','diagQuadratic','prior',[0.5 0.5])));
              Xm1=reorder_D_cols(Xm,find(~ismember(Xm.supdat(end,:),tst)));
              Xm1.supdat(end,:)=Xm1.supdat(end,:)-1;
              %      Xm1=Dtrain; % needed to select from all genes
              global FET_SEL_SI
              FET_SEL_SI=cell(size(Xm1.dat,2),D_n_cls(Xm1,size(Xm1.supdat,1)));
              [confmat,nerr,cvres,D,params]=crossvalidate(Xm1,...
                                                          size(Xm1.supdat,1),...
                                                          struct('method','loocv','classifier',...
                                                                struct('method','gg_pnn_w_sel_fet',...
                                                                'verbose',1,...
                                                                'num_features',size(Xm.dat,1),...
                                                                'fet_sel',struct('method','gp',...
                                                                'score',struct('method','ttest_minvar','minvar',0.3^2)),...
                                                                'sig_type',2,...
                                                                'sig',-1,...
                                                                'dist_type','cosine','class_weights',[0.5 0.5]')));
            end

        end
        
        % gg_pnn: crossprod_struct(struct('num_features',
        % num2cell(10:10:300)),struct('sig',num2cell([-1:-0.5:-4])))
        %         cosine, euclid
        %         gp  2class_pv
        % nb: struct('num_features', num2cell(10:10:300)) gp 2class_pv
        
        params=struct('method','cv_select',...
                      'hyper',crossprod_struct(struct('num_features', num2cell(10:10:300)),struct('sig',num2cell([-1:-0.5:-4]))),...
                      'show_performance',0,...
                      'cv',struct('method','loocv',...
                                  'classifier',struct('method','classifier_w_sel_fet',...
                                                      'classifier_type','gg_pnn',...
                                                      'fet_sel',struct('method','gp',... %2class_pv',... %'gp' '2class'
                                                          'score',struct('method','ttest_minvar', ...
                                                          'minvar',(0.3)^2)),...
                                                      'sig_type',2,'sig',-1,'class_weights',[0.5; 0.5],...
                                                      'dist_type','cosine')),...
                      'classifier',struct('method','classifier_w_sel_fet',...
                                          'classifier_type','gg_pnn',...
                                          'fet_sel',struct('method','gp',... %'2class_pv',... %'gp' '2class'
                                                          'score',struct('method','ttest_minvar', ...
                                                          'minvar',(0.3)^2)),...
                                          'sig_type',2,'sig',-1,'class_weights',[0.5; 0.5],...
                                          'dist_type','cosine'));
        for t=t_vals
          switch t
           case 1
            params.hyper=crossprod_struct(struct('num_features', num2cell(10:10:300)),struct('sig',num2cell([-1:-0.5:-4]))); %4
            params.cv.classifier.classifier_type='gg_pnn';
            params.cv.classifier.fet_sel.method='gp';
            params.cv.classifier.dist_type='cosine';
            params.classifier.classifier_type='gg_pnn';
            params.classifier.fet_sel.method='gp';
            params.cv.classifier.dist_type='cosine';
            
           case 2
            params.hyper=crossprod_struct(struct('num_features', num2cell(10:10:300)),struct('sig',num2cell([-1:-0.5:-4]))); %4
            params.cv.classifier.classifier_type='gg_pnn';
            params.cv.classifier.fet_sel.method='gp';
            params.cv.classifier.dist_type='euclidean';
            params.classifier.classifier_type='gg_pnn';
            params.classifier.fet_sel.method='gp';
            params.cv.classifier.dist_type='euclidean';
            
           case 3
            params.hyper=crossprod_struct(struct('num_features', num2cell(10:10:300)),struct('sig',num2cell([-1:-0.5:-4]))); %4
            params.cv.classifier.classifier_type='gg_pnn';
            params.cv.classifier.fet_sel.method='2class_pv';
            params.cv.classifier.dist_type='cosine';
            params.classifier.classifier_type='gg_pnn';
            params.classifier.fet_sel.method='2class_pv';
            params.cv.classifier.dist_type='cosine';
            
           case 4
            params.hyper=crossprod_struct(struct('num_features', num2cell(10:10:300)),struct('sig',num2cell([-1:-0.5:-4]))); %4
            params.cv.classifier.classifier_type='gg_pnn';
            params.cv.classifier.fet_sel.method='2class_pv';
            params.cv.classifier.dist_type='euclidean';
            params.cv.keep_cvres=1;
            params.classifier.classifier_type='gg_pnn';
            params.classifier.fet_sel.method='2class_pv';
            
           case 5
            params.hyper=struct('num_features', num2cell(10:10:300));
            params.cv.classifier.classifier_type='matlab_classify';
            params.cv.classifier.prior=[0.5; 0.5];
            params.cv.classifier.dist_type='diagQuadratic';          
            params.cv.classifier.fet_sel.method='gp';
            params.classifier.classifier_type='matlab_classify';
            params.classifier.fet_sel.method='gp';
            params.classifier.prior=[0.5; 0.5];
            params.classifier.dist_type='diagQuadratic';          
            
           case 6
            params.hyper=struct('num_features', num2cell(10:10:300)); % 300
            params.cv.classifier.classifier_type='matlab_classify';
            params.cv.classifier.prior=[0.5; 0.5];
            params.cv.classifier.dist_type='diagQuadratic';          
            params.cv.classifier.fet_sel.method='2class_pv';
            params.cv.keep_cvres=1;
            params.classifier.classifier_type='matlab_classify';
            params.classifier.fet_sel.method='2class_pv';
            params.classifier.prior=[0.5; 0.5];
            params.classifier.dist_type='diagQuadratic';          
           
           case 7
            params.hyper=struct('sig',num2cell([-1:-0.5:-4])); %4
            params.cv.classifier.classifier_type='gg_pnn_pca';
            params.cv.classifier.pca_k=3;  
            params.cv.classifier.num_features=-1;              
            params.cv.classifier.fet_sel.method='2class_pv_fdr';
            params.cv.classifier.fet_sel.fdr_thresh=0.1;
            params.cv.classifier.dist_type='euclidean';
            params.classifier.classifier_type='gg_pnn_pca';
            params.classifier.num_features=-1;              
            params.classifier.pca_k=3;            
            params.classifier.fet_sel.method='2class_pv_fdr';
            params.classifier.fet_sel.fdr_thresh=0.1;
            params.classifier.dist_type='euclidean';
                       
           case 8
            params.hyper=struct('num_features', num2cell(10:10:300)); % 300
            params.cv.classifier.classifier_type='matlab_classify';
            params.cv.classifier.prior=[0.5; 0.5];
            params.cv.classifier.dist_type='diaglinear';          
            params.cv.classifier.minvar=(0.3)^2;
            params.cv.classifier.fet_sel.method='2class_pv';
            params.cv.keep_cvres=1;
            params.classifier.classifier_type='matlab_classify';
            params.classifier.fet_sel.method='2class_pv';
            params.classifier.prior=[0.5; 0.5];
            params.classifier.dist_type='diaglinear';          
            params.classifier.minvar=(0.3)^2;
             
           case 9
            params.hyper=struct('num_features', num2cell(10:10:300)); % 300
            params.cv.classifier.classifier_type='matlab_classify';
            params.cv.classifier.prior=[0.5; 0.5];
            params.cv.classifier.dist_type='diaglinear_median';          
            params.cv.classifier.minvar=(0.3)^2;
            params.cv.classifier.fet_sel.method='2class_pv';
            params.cv.keep_cvres=1;
            params.classifier.classifier_type='matlab_classify';
            params.classifier.fet_sel.method='2class_pv';
            params.classifier.prior=[0.5; 0.5];
            params.classifier.dist_type='diaglinear_median';          
            params.classifier.minvar=(0.3)^2;
             
          end 
          classify_and_cluster(X,Dtrain,Dtest,phens,cur_supid,params,t,lesion_name,options,i, ...
                                               force_run);
        end
        
        if (0)
          %    params.xlsfdir='/xchip/data/gadgetz/lsfres/';
          %    params.lsf_maxjobs=11;
          params=struct('method','cv_select',...
                          'hyper',crossprod_struct(struct('num_features', num2cell(10:10:300)),struct('sig',num2cell([-1:-0.5:-4]))),...
                          'show_performance',0,...
                          'cv',struct('method','loocv',...
                                      'classifier',struct('method','classifier_w_sel_fet',...
                                                          'classifier_type','gg_pnn',...
                                                          'fet_sel',struct('method','gp',... %2class_pv',... %'gp' '2class'
                                                                           'score',struct('method','ttest_minvar', ...
                                                                                          'minvar',(0.3)^2)),...
                                                          'sig_type',2,'sig',-1,'class_weights',[0.5; 0.5],...
                                                          'dist_type','cosine')),...
                          'classifier',struct('method','classifier_w_sel_fet',...
                                              'classifier_type','gg_pnn',...
                                              'fet_sel',struct('method','gp',... %'2class_pv',... %'gp' '2class'
                                                               'score',struct('method','ttest_minvar', ...
                                                                              'minvar',(0.3)^2)),...
                                              'sig_type',2,'sig',-1,'class_weights',[0.5; 0.5],...
                                              'dist_type','cosine'));

                                          
            params=struct('method','cv_select',...
                          'hyper',struct('num_features', num2cell(10:10:300)),...
                          'show_performance',0,...
                          'cv',struct('method','loocv',...
                                      'classifier',struct('method','classifier_w_sel_fet',...
                                                          'classifier_type','matlab_classify',...
                                                          'fet_sel',struct('method','2class_pv',... %'gp' '2class'
                                                                           'score',struct('method','ttest_minvar', ...
                                                                                          'minvar',(0.3)^2)),...
                                                          'prior',[0.5; 0.5],...
                                                          'dist_type','diagQuadratic')),...  %'diagQuadratic'
                          'classifier',struct('method','classifier_w_sel_fet',...
                                              'classifier_type','matlab_classify',...
                                              'fet_sel',struct('method','2class_pv',... %'2class'
                                                               'score',struct('method','ttest_minvar', ...
                                                                              'minvar',(0.3)^2)),...
                                              'prior',[0.5; 0.5],...
                                              'dist_type','diagQuadratic'));
                                                           

            global FET_SEL_SI
            FET_SEL_SI=[] % cell(size(Dtrain.dat,2)); %,D_n_cls(Dtrain,size(Dtrain.supdat,1)));
            [pred,tr,res,cls,Dtest1,Dtrain1]=train_and_predict(Dtrain,Dtest,size(Dtrain.supdat,1),params);
            
            Xm=reorder_D_sup(reorder_D_rows(X,cls.cls.si),'cols',[ find_supid(X,phens) cur_supid]);

            [Xm_ord,Xm_ord.gdend,Xm_ord.sdend]=two_way_clustering(Xm,...
                                                              struct('dist','correlation','preproc',...
                                                              'row_center_and_normalize','cluster','average'),...
                                                              struct('dist','euclidean','preproc',...
                                                              'row_center_and_normalize','cluster','ward'));
            
            
            display_D(Xm_ord,Xm_ord.gdend,Xm_ord.sdend,'dends');
            print_D([ lesion_name '_' options{i,1} '.Naive_Bayes.xv.clustering'],{{'pdf'}},1);
            disp([ lesion_name ' ' options{i,1} ' Naive_Bayes xv: ' num2str(length(m1)) ' genes']);
            
            
        end
    end
end


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GSEA

XS=collapse_to_symbols(X,'median');

exp_fname=['gene_expression.symbol.no' chrst '.' lesion_name '.' num2str(size(X.dat,2)) '.gct'];
if ~exist(exp_fname,'file')
    write_mit_gct_file(exp_fname,XS);
end
write_mit_cls_file([lesion_name '.' num2str(size(X.dat,2)) '.cls'],X,find_supid(X,lesion_name),1,1);

XSf=filter_D_rows(XS,struct('method','std','thresh',0.3));

gene_set_files={'~/gsea/rameen_sets/final_sets/rameen.symbol.gmt',...
    '~/gsea/sets_060210/c2.symbols.gmt',...
    '~/gsea/sets_060210/xhx_and_bartel.symbols.gmt'};

for i=1:length(gene_set_files)
    sr{i}=read_mit_gmt_file(gene_set_files{i});
    XSfs{i}=gsea_match_sets(sr{i},XSf,1);
end

cd([ glad_dir '/with_GE/' alg_dir '/' subset_dir '/' dir_name '/' lesion_name]);


% Gx=add_D_sup(G,C.supacc,C.supdesc,C.supdat,'cols');


for i=1:length(gene_set_files)
    XP{i}=project_D_on_gsupdat(XSfs{i},'sum_of_normalized');
end

G=reorder_D_sup(G,'cols',find_supid(G,phens));
[dum,si]=sortrows(G.supdat');

G=reorder_D_cols(G,si);

for i=1:length(gene_set_files)
    YPp{i}=reorder_D_cols(reorder_D_sup(YP{i},'cols',find_supid(YP{i},phens,'exact')),si);
end
Cp=reorder_D_cols(C,si);
clear C
clear G Gx

if two_class
    options={'No_versus_Yes',1,2};
else
    options={'No_versus_High',1,3;'No_versus_Low',1,2;'Low_versus_High',2,3;'No_versus_REST',1,[2 3]};
end
% options={'No_versus_REST',1,[2 3]};

for i=1:size(options,1)
    opt_dir=[ glad_dir '/with_GE/' alg_dir '/' subset_dir '/' dir_name '/' lesion_name '/' options{i,1}];
    disp(opt_dir);
    if ~exist(opt_dir,'dir')
        mkdir(opt_dir);
    end
    cd(opt_dir);
    use_idx=find(ismember(YPp{1}.supdat(find_supid(YPp{1},lesion_name),:),[ options{i,2} options{i,3}]));
    GpU=reorder_D_cols(Gp,use_idx);
    for j=1:length(gene_set_files)
        YPpU{j}=reorder_D_cols(YPp{j},use_idx);
    end
    CpU=reorder_D_cols(Cp,use_idx);

    t=0.3;
    calls=call_regs(CpU,regs,t);
    calls1={ismember(YPpU{1}.supdat(find_supid(YPp{1},lesion_name),:),options{i,3})};
    if isempty(regid)
        disp('Using REG #1 ...');
        regid=1;
    end
    regs1={regs{ad}(regid)};
    pvs1={pvs{ad}};

    regs1{1}
    if std(calls1{1},0,2)>0
        tab=snp_generate_region_data(CpU,GpU,...
            strvcat(Ys{1}.ll),rg,regs1,calls1,t,pvs1,...
            YPpU);
        clear CpU
        %bsub -J job1 /usr/java/java1.4/bin/java -cp /xchip/data/gadgetz/gsea_jar/gsea2.jar -Xmx1800m xtools.gsea.Gsea
        %-param_file ~/gsea/gsea.weighted.prm -metric tTest -rpt_label No_vs_High_tTest -gmx ~/gsea/rameen_sets/final_sets/rameen.symbol.gmt -min_set 3
        %-res gene_expression.symbol.no7.62.gct -cls EGFR_amp.62.cls#No_versus_High -out ~/projects/snp/Rameen/new/glad_dChip/with_GE/ALL/gsea

        cur_dir=pwd;
        for j=1:2 % length(gene_set_files)
            unix_str=[ 'bsub -J gsea /usr/java/jdk1.5.0_03/jre/bin/java -cp /xchip/projects/xtools/xtools.jar -Ddebug=true ' ...
                '-Xmx1800m xtools.gsea.Gsea ' ...
                'gsea -param_file ~/gsea/gsea.weighted.prm -metric Signal2Noise ' ...
                '-rpt_label ' options{i,1} '_S2N_' regexprep(regexprep(gene_set_files{j},'.*/',''),'\..*','') ...
                ' -gmx ' gene_set_files{j} ' ' ...
                '-res ../gene_expression.symbol.no' chrst '.' lesion_name '.' num2str(size(CpU.dat,2)) '.gct -cls ../' lesion_name '.' num2str(size(CpU.dat,2)) '.cls#' ...
                options{i,1} ' -set_min 3 ' ...
                '-out ' pwd ];
            if (0)
                unix_str=[ '/xchip/projects/xtools/run bsub gsea -param_file ~/gsea/gsea.weighted.prm -metric Signal2Noise ' ...
                    '-rpt_label ' options{i,1} '_S2N_' regexprep(regexprep(gene_set_files{j},'.*/',''),'\..*','') ...
                    ' -gmx ' gene_set_files{j} ' ' ...
                    '-res ../gene_expression.symbol.no' chrst '.' lesion_name '.62.gct -cls ../' lesion_name '.62.cls#' options{i,1} ' -set_min 3 ' ...
                    '-out ' pwd ];
            end
            disp(unix_str);
            if ~exist('no_gsea','var') || ( exist('no_gsea','var') && ~no_gsea)
                unix(unix_str);
            end
        end
    end
    clear CpU
end


%------------------------------------------------------------------------------------------------------

function classify_and_cluster(X,Dtrain,Dtest,phens,cur_supid,params,t,lesion_name,options,i,force_run)
global FET_SEL_SI
FET_SEL_SI=[] % cell(size(Dtrain.dat,2)); %,D_n_cls(Dtrain,size(Dtrain.supdat,1)));

fdr_mat_fname=[ lesion_name '_' options{i,1} '.' num2str(t) '.fdr_res.mat'];
if force_run || ~exist(fdr_mat_fname,'file')
  fdr_level=0.05;
  sel_method=struct('method','bhfdr','thresh',fdr_level);
  [M.idx,M.q,M.p,M.s,M.pi0,M.F]=get_top_markers(Dtrain,cur_supid,params.classifier.fet_sel.score,...
                                                       -1,...
                                                       sel_method);
  [sq,si]=sort(M.q);
  M.idx=M.idx(si);
  M.q=M.q(si);
  M.p=M.p(si);
  M.s=M.s(si);
  M.score=cur_supid,params.classifier.fet_sel.score;
  M.sel_method=sel_method;
  save(fdr_mat_fname,'M');
else
  load(fdr_mat_fname);
end

X=add_D_sup(X,['FDR_' num2str(M.sel_method.thresh)],['FDR_' num2str(M.sel_method.thresh)],ones_in(size(X.dat,1), ...
                                                  M.idx),'rows');
X.gscore=M.F.s;
X.gq=M.F.q;
X.gp=M.F.p;
X=add_D_field(X,'gene',{'gscore','gq','gp'});


if ~isempty(M.idx)
  T=reorder_D_rows(Dtrain,M.idx);
  T.gscore=M.s;
  T.gq=M.q;
  T.gp=M.p;
  write_D_genes([ lesion_name '_' options{i,1} '.' num2str(t) '.fdr.' num2str(M.sel_method.thresh) '.' ...
                  num2str(length(M.idx)) '.genes.txt'],...
                T,{'gscore','gq','gp'});
end

pred_mat_fname=[ lesion_name '_' options{i,1} '.' num2str(t) '.pred_res.mat'];

if force_run || ~exist(pred_mat_fname,'file')
  [pred,tr,res,cls,Dtest1,Dtrain1]=train_and_predict(Dtrain,Dtest,size(Dtrain.supdat,1),params);
  save(pred_mat_fname,'pred','tr','cls');
else
  load(pred_mat_fname);
end
disp(['Case ' num2str(t) ' Res = ' num2str(pred')]);


sup_name=[ lesion_name '_' options{i,1} '.' num2str(t) ];
X=add_D_sup(X,sup_name,sup_name,ones_in(size(X.dat,1),cls.cls.si),'rows');

report_gene_sets([ lesion_name '_' options{i,1} '.' num2str(t) '.gene_sets.txt'],X,1:(size(X.gsupdat,1)-2), ...
                 (size(X.gsupdat,1)-1):size(X.gsupdat,1));


Xm=reorder_D_sup(reorder_D_rows(X,cls.cls.si),'cols',[ find_supid(X,phens) cur_supid]);
disp([ lesion_name ' ' options{i,1} ' ' num2str(t) ' xv: ' num2str(length(cls.cls.si)) ' genes']);

for cvi=1:length(cls.f_res)
  nerr(cvi)=cls.f_res(cvi).confmat(1,2)+cls.f_res(cvi).confmat(2,1);
  fp(cvi)=cls.f_res(cvi).confmat(2,1);
  mls(cvi)=cls.f_res(cvi).mean_logscore;
  if isfield(cls.params.hyper(cvi),'num_features')
    nf(cvi)=cls.params.hyper(cvi).num_features;
  else
    nf(cvi)=NaN; % put number of FDR!
  end
end
maxnerr=max(nerr);
nsamples=sum(sum(cls.f_res(1).confmat));
mls(isinf(mls))=NaN;
%mls=mls./nf;

fs=14;
figure(1); clf
set(gcf,'Visible','off');
plot(nf,nerr,'o-'); hold on
plot(nf,fp,'o--');
plot(nf,(mls-min(mls))./(max(mls)-min(mls))*nsamples,'.-k');
ax=axis;
axis([ax(1:2) 0 nsamples]);
minerr=find(nerr==min(nerr));
[tmp,choice]=max(mls(minerr));
choice=minerr(choice);
ax=axis;
line([ nf(choice) nf(choice)],ax(3:4),'Color','r');
xlabel('Number of features','FontSize',fs);
ylabel('Number of errors','FontSize',fs);
legend({'N_{err}','FP','MLS'});
print_D([ lesion_name '_' options{i,1} '.' num2str(t) '.xv.params'],{{'pdf'}},1);


[Mtrn,mtrn1,mtrn2]=match_string_sets(Dtrain.sdesc,X.sdesc);
[Mtst,mtst1,mtst2]=match_string_sets(Dtest.sdesc,X.sdesc);

v=nan(1,size(X.dat,2));
v(mtrn2)=Dtrain.supdat(cur_supid,mtrn1);
v(mtst2)=Dtest.supdat(cur_supid,mtst1);
v(v==options{i,4})=nan;

if ismember(t,[6 8 9])
  % from NB
  w=nan(1,size(X.dat,2));
  w(mtst2)=cls.cls.posterior(mtst1,2)';
  for j=1:length(mtrn1)
    w(mtrn2(j))=cls.f_res(choice).cvres(mtrn1(j)).classifier.posterior(2);
  end
else % for PNN
  w=nan(1,size(X.dat,2));
  w(mtst2)=cls.cls.score(mtst1);
  for j=1:length(mtst2)
    if ~pred(j)
      w(mtst2(j))=1-w(mtst2(j));
    end    
  end  
  keyboard
  for j=1:length(mtrn1)
    w(mtrn2(j))=cls.f_res(choice).cvres(mtrn1(j)).classifier.score;
    if ~cls.f_res(choice).cvres(mtrn1(j)).pred
      w(mtrn2(j))=1-w(mtrn2(j));
    end
  end
end

x=options{i,2}*ones(1,size(X.dat,2));
x(w>0.5)=options{i,3};
x(mtrn2)=x(mtrn2).*(-1).^(x(mtrn2)~=v(mtrn2));

comp_str=[lesion_name '_' options{i,1}];
sup_str=[repmat(comp_str,3,1) strvcat('_TRN','_POST','_PRED')];
[X,add_supid]=add_D_sup(X,sup_str,sup_str,[v;w;x],'cols');
Xtmp=reorder_D_sup(X,'cols',[find_supid(X,phens) add_supid]);
write_eisen_dat([comp_str '.' num2str(t) '.prediction.txt' ],strvcat(Xtmp.supacc),strvcat(Xtmp.supdesc),Xtmp.sdesc, ...
                'PREDICTION',Xtmp.supdat,[],[],0);

clustering_m={'ward','average','centroid'};
for j=1:length(clustering_m)
  [Xm_ord,Xm_ord.gdend,Xm_ord.sdend]=two_way_clustering(Xm,...
                                                    struct('dist','correlation','preproc',...
                                                    'row_center_and_normalize','cluster','average'),...
                                                    struct('dist','euclidean','preproc',...
                                                    'row_center_and_normalize','cluster',clustering_m{j}));
  
  Xm_ord.gacc=strcat(Xm_ord.gacc,cellstr(repmat(':',length(Xm_ord.gacc),1))',Xm_ord.gsymb);
  figure(1); clf
  set(gcf,'Visible','off');
  Xm_ord=fix_supdat(Xm_ord);
  display_D(Xm_ord,Xm_ord.gdend,Xm_ord.sdend,'dends');
  print_D([ lesion_name '_' options{i,1} '.' num2str(t) '.xv.' clustering_m{j} '.clustering'],{{'pdf'}},1);
  disp(clustering_m{j});
end
write_D_genes([ lesion_name '_' options{i,1} '.' num2str(t) '.' num2str(size(Xm_ord.dat,1)) '.genes.txt'],Xm_ord,{'gscore','gq','gp'});








