function [res,fullres]=gp_pnn(train_gct,train_cls,test_gct,test_cls,fullres_fname,predres_fname,num_fets,...
    test_type,min_std,sig_type,sig,prior_type,prior,dist_type,load_params,params_file)
% gp_pnn(train_gct,train_cls,test_gct,test_cls,predres,sig,sig_type,prior_ratio,dist_type)
% 

% read train
if ischar(train_gct)
    Dtrain=read_mit_gct_file(train_gct);
    Dtrain=read_mit_cls_file(Dtrain,train_cls,1);
else
    Dtrain=train_gct;
end
n_cls=D_n_cls(Dtrain,1);
n_train=size(Dtrain.dat,2);

% read test
if ischar(test_gct)
    Dtest=read_mit_gct_file(test_gct);
    Dtest=read_mit_cls_file(Dtest,test_cls,1);
else
    Dtest=test_gct;
end

n_test=size(Dtest.dat,2);

% handle inputs
if ischar(load_params)
    load_params=str2num(load_params);
end

if load_params
    params_odf=read_mit_odf_file(params_file);
    num_fets=params_odf.data{2};
    test_type=params_odf.data{3};
    min_std=params_odf.data{4};
    sig_type=params_odf.data{5};
    sig=params_odf.data{6};
    prior_type=params_odf.data{7};
    prior=params_odf.data{8};
    dist_type=params_odf.data{9};
else
    if ischar(num_fets)
        num_fets=str2num(num_fets);
    end
    
    if ischar(test_type)
        test_type=str2num(test_type);
    end
    
    if ischar(min_std)
        min_std=str2num(min_std);
    end
    
    if ischar(sig_type)
        sig_type=str2num(sig_type);
    end
    
    if ischar(sig)
        sig=str2num(sig);
    end
    
    if ischar(prior_type)
        prior_ratio=str2num(prior_type);
    end  
    
    if ischar(prior)
        prior=str2num(prior);
    end
    
    dist_type_num=str2num(dist_type);
    if ~isempty(dist_type_num)
        switch(dist_type_num)
            case 0
                dist_type='euclidean';
            case 1
                dist_type='cosine';
            otherwise
                error('No such distance type.');
        end
    end
    
    num_fets=repmat(num_fets,n_cls,1);
    test_type=repmat(test_type,n_cls,1);
    min_std=repmat(min_std,n_cls,1);
    sig_type=repmat(sig_type,n_cls,1);
    sig=repmat(sig,n_cls,1);
    prior_type=repmat(prior_type,n_cls,1);
    prior=repmat(prior,n_cls,1);
    dist_type=cellstr(repmat(dist_type,n_cls,1));
end

used_fets={};
scores=[];
logscores=[];
used_sigs=[];
pred_classes=[];
true_classes=[];
correct=[];
cvec=[];
for c=1:n_cls
    train_cls=[ double(Dtrain.supdat(1,:)~=c); double(Dtrain.supdat(1,:)==c) ]; 
    test_cls=[ double(Dtest.supdat(1,:)~=c); double(Dtest.supdat(1,:)==c) ]; 
    
    % Select features
    switch(test_type(c))
        case 0
            diff_test_type='gcsnr';
        case 1
            diff_test_type='ttest';
        case 2
            diff_test_type=struct('method','ttest_minvar','minvar',min_std(c).^2);
        otherwise
            error('Test type not supported.');
    end
    
    global FET_SEL_SI
    if ~isfield(Dtrain,'fet_sel_id') || isempty(FET_SEL_SI{Dtrain.fet_sel_id,c})
        [dum,s]=differential_analysis(Dtrain,find(train_cls(1,:)),find(train_cls(2,:)),diff_test_type,1);
        [ss,si]=sort(s);
        if isfield(Dtrain,'fet_sel_id')
            FET_SEL_SI{Dtrain.fet_sel_id,c}=si;
%            fprintf(1,'x');
        end
    else
        si=FET_SEL_SI{Dtrain.fet_sel_id,c};
%        fprintf(1,'.');
    end
    
    nfet_top=ceil(num_fets(c)/2);
    nfet_bot=num_fets(c)-nfet_top;
    si_top=si((end-nfet_top+1):end);
    si_bot=si(1:nfet_bot);
    use_fets=flipud([si_bot; si_top]);
    
    Dtrain_c=reorder_D_rows(Dtrain,use_fets);
    Dtest_c=reorder_D_rows(Dtest,use_fets);
    
    % Set prior
    if prior_type(c)==0 % empirical, based on train set
        prior(c)=sum(train_cls(2,:),2)/n_train;
    end
    
    w=[ 1-prior(c); prior(c)]; 
    
    % classify
    [class,score,used_sig,logscore]=pnn_classify_standalone(Dtest_c.dat',Dtrain_c.dat',...
        train_cls,...
        w,sig_type(c),sig(c),dist_type{c},...
        test_cls);
    class=class-1;
    st=sprintf('%s,',Dtest_c.gacc{:});
    used_fets(end+1)={st(1:(end-1))};
    if n_test>1
        used_fets((end+1):(end+n_test-1))=cellstr(repmat('.',n_test-1,1));
    end
    scores=[scores score];
    logscores=[logscores logscore];
    used_sigs=[used_sigs  used_sig];
    pred_classes=[pred_classes class];
    true_classes=[true_classes test_cls(2,:)];
    correct=[correct class==test_cls(2,:)];
    cvec=[cvec c*ones(1,length(class))];
end

cls1_scores=scores;
cls1_scores(pred_classes==0)=1-cls1_scores(pred_classes==0);
scoremat=reshape(cls1_scores,n_test,n_cls)';
[maxscore,pred_class]=max(scoremat,[],1);

% write results

res.col_names={ 'Sample' 'True class' 'Predicted class' 'Score' 'Correct?'};
res.col_types={ 'string' 'string' 'string' 'float' 'boolean'};
res.data={ cellstr(Dtest.sdesc) ...
        cellstr(num2str(Dtest.supdat(1,:)')) cellstr(num2str(pred_class')) ...
        maxscore (Dtest.supdat(1,:)==pred_class)'};
res.Model='Prediction Results';
res.PredictorModel='PNN';
res.NumCorrect=num2str(sum(Dtest.supdat(1,:)==pred_class));
res.NumErrors=num2str(n_test-str2num(res.NumCorrect));

if ~isempty(predres_fname)
  write_mit_odf_file(predres_fname,res);
end

fullres.used_sigmas=used_sigs;
fullres.col_names={ 'Sample' 'True class' 'Prediction for class' ...
        'True in-class' 'Predicted in-class' 'Score' 'Logscore' 'Correct?' 'Used features'};
fullres.col_types={ 'string' 'string' 'string' 'string' 'string' 'float' 'float' 'boolean' 'string'};
fullres.data={ repmat(cellstr(Dtest.sdesc),n_cls,1) ...
        repmat(cellstr(num2str(Dtest.supdat(1,:)')),n_cls,1) cellstr(num2str(cvec')) ...
        cellstr(num2str(true_classes')) cellstr(num2str(pred_classes')) scores logscores correct used_fets};
fullres.Model='PNN';
if ~isempty(fullres_fname)
    write_mit_odf_file(fullres_fname,fullres);
end

% mcc -m gp_pnn





