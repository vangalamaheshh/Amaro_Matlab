function [pred,tr,res,cls,Dtest,Dtrain]=train_and_predict(Dtrain,Dtest,supid,classifier_params)

if ischar(classifier_params)
  classifier_params.method=classifier_params;
end

cls.params=classifier_params;
switch classifier_params.method
 case 'cv_select'
  fres=1:length(classifier_params.hyper);
  classifier_params.cv.presend_files=1;
  global FET_SEL_SI
  FET_SEL_SI=cell(size(Dtrain.dat,2),1);
  for f=1:size(fres,2)
    if isfield(classifier_params.cv,'params')
      classifier_params.cv.params= ...
          add_struct(classifier_params.cv.params,classifier_params.hyper(f));
    else
      classifier_params.cv.params=classifier_params.hyper(f);
    end
    classifier_params.cv.params
    [f_confmat,f_nerr,f_cvres,Dtrain,params]=crossvalidate(Dtrain, ...
                                                      supid, ...
                                                      classifier_params.cv);
    if ~isempty(params)
      classifier_params.cv.presend_files=0;
      classifier_params.cv.params= ...
          add_struct(classifier_params.cv.params,params);
    end      
    fres(2,f)=f_nerr;
    f_res(f).confmat=f_confmat;
    if isfield(classifier_params.cv,'keep_cvres') && ...
          classifier_params.cv.keep_cvres
      
      f_res(f).cvres=f_cvres;
    end
    tmp=cat(1,f_cvres(:).classifier);
    f_res(f).mean_logscore=mean(cat(1,tmp(:).logscore));
    fres(3,f)=f_res(f).mean_logscore;
    Dtrain=append_D_sup(Dtrain,size(Dtrain.supdat,1),num2str(fres(1,f)),...
                        [' ' num2str(fres(1,f))]);
    disp([fres(1,f) f_nerr f_res(f).mean_logscore ]);
  end
  if exist('ne_lms.mat')
    zzz=load('ne_lms.mat');
    zne=zzz.zne;
    zmls=zzz.zmls;
    zne=[zne; fres(2,:)];
    zmls=[zmls; fres(3,:)];
    save ne_lms.mat zne zmls
  else
    zne=[fres(2,:)];
    zmls=[fres(3,:)];
    save ne_lms.mat zne zmls
  end
  [min_err,f_choose_i]=min(fres(2,:));
  with_min_err=find(fres(2,:)==min_err);
  if length(with_min_err)>1
    [max_logscore,maxi]=max(fres(3,with_min_err));
    f_choose_i=with_min_err(maxi);
  end
  FET_SEL_SI=[];
  cls.f_vals_chosen=classifier_params.hyper(f_choose_i);
  cls.fres=fres;
  cls.f_res=f_res;
  cls.min_err=min_err;
  Dtrain=append_D_sup(Dtrain,size(Dtrain.supdat,1)-size(fres,2)+f_choose_i,'*',' *');
    
  if isfield(classifier_params,'show_performance') && classifier_params.show_performance
    figure(1); clf
    plot(fres(1,:),fres(2,:),'x');
    drawnow;
%    disp('hit any key');
%    pause;
  end
  if isfield(classifier_params.classifier,'params')
    classifier_params.classifier.params= ...
        add_struct(classifier_params.cv.params,cls.f_vals_chosen);
  else
    classifier_params.classifier.params=cls.f_vals_chosen;
  end
  classifier_params.classifier.params
  [pred,tr,res,cls.cls,Dtest,Dtrain]=train_and_predict(Dtrain,Dtest,supid,classifier_params.classifier);

  
 case 'cv_select_features'
  if isfield(classifier_params,'f_vals')
    fres=classifier_params.f_vals;
  else
    fres=classifier_params.f_start:classifier_params.f_step: ...
         min(classifier_params.f_end,size(Dtrain.dat,1));
  end
  classifier_params.cv.presend_files=1;
  for f=1:size(fres,2)
    classifier_params.cv.params.num_features=num2str(fres(1,f));
    [f_confmat,f_nerr,f_cvres,Dtrain,params]=crossvalidate(Dtrain, ...
                                                      supid, ...
                                                      classifier_params.cv);
    if ~isempty(params)
      classifier_params.cv.presend_files=0;
      classifier_params.cv.params= ...
          add_struct(classifier_params.cv.params,params);
    end      
    fres(2,f)=f_nerr;
    f_res(f).confmat=f_confmat;
    f_res(f).cvres=f_cvres;
    tmp=cat(1,f_cvres(:).classifier);
    f_res(f).mean_logscore=mean(cat(1,tmp(:).logscore));
    fres(3,f)=f_res(f).mean_logscore;
    Dtrain=append_D_sup(Dtrain,size(Dtrain.supdat,1),num2str(fres(1,f)),...
                        [' ' num2str(fres(1,f))]);
    disp([fres(1,f) f_nerr f_res(f).mean_logscore ]);
  end
  [min_err,f_choose_i]=min(fres(2,:));
  with_min_err=find(fres(2,:)==min_err);
  if length(with_min_err)>1
    [max_logscore,maxi]=max(fres(3,with_min_err));
    f_choose_i=with_min_err(maxi);
  end
  cls.f_n_chosen=fres(1,f_choose_i);
  cls.fres=fres;
  cls.f_res=f_res;
  cls.min_err=min_err;
  Dtrain=append_D_sup(Dtrain,size(Dtrain.supdat,1)-size(fres,2)+f_choose_i,'*',' *');
    
  if isfield(classifier_params,'show_performance') && classifier_params.show_performance
    figure(1); clf
    plot(fres(1,:),fres(2,:),'x');
    drawnow;
%    disp('hit any key');
%    pause;
  end
  classifier_params.classifier.params.num_features= ...
      num2str(cls.f_n_chosen);
  [pred,tr,res,cls.cls,Dtest,Dtrain]=train_and_predict(Dtrain,Dtest,supid,classifier_params.classifier);

 case 'one_vs_rest'
  tr=Dtest.supdat(supid,:)';
  class_types=unique(Dtrain.supdat(supid,:));
  [typeacc,typedesc,Dtrain,range,non_empty]=decollapse_supdat(Dtrain,1);
  [typeacc_tr,typedesc_tr,Dtest,range_tr,non_empty_tr]=decollapse_supdat(Dtest,1,non_empty);
  cls.class_types=non_empty;
  cls.range=range;
  cls.pred=[];  
  
  if isfield(classifier_params,'lsfdir')
    l=lsf(classifier_params.lsfdir);
    l=maxjobs(l,classifier_params.lsf_maxjobs);
    h=zeros(length(class_types),1);
    for i=1:length(class_types)
      res(i).class_type=class_types(i);
      [l,h(i)]=bsub(l,{'pred','tr','res','cls','Dtest','Dtrain'},'train_and_predict',...
                    {Dtrain,Dtest,range(find(non_empty==class_types(i))),classifier_params.classifier});
    end
    [l,lres]=wait(l);
    for i=1:length(class_types)
      res(i).pred=lres{h(i)}.pred;
      res(i).tr=lres{h(i)}.tr;
      res(i).res=lres{h(i)}.res;
      res(i).cls=lres{h(i)}.cls;
      res(i).Dtest=lres{h(i)}.Dtest;
      res(i).Dtrain=lres{h(i)}.Dtrain;
      cls.pred=[cls.pred; res(i).pred'];
    end
  else
    for i=1:length(class_types)
      res(i).class_type=class_types(i);
      [res(i).pred,res(i).tr,res(i).res,res(i).cls,res(i).Dtest,res(i).Dtrain]=train_and_predict(Dtrain,Dtest, ...
                                                        range(find(non_empty==class_types(i))),classifier_params.classifier);
      cls.pred=[cls.pred; res(i).pred'];
    end
  end
  for i=1:length(res)
    cls.score_mat(i,:)=res(i).pred'.*res(i).cls.cls.score + (1-res(i).pred').*(1-res(i).cls.cls.score);
  end
  [cls.score,pred]=max(cls.score_mat,[],1);
  pred=cls.class_types(pred');
  if ~isfield(classifier_params,'keep_data') || ...
        ~classifier_params.keep_data
    if isfield(res,'Dtest') 
      res=rmfield(res,'Dtest');
    end
    if isfield(res,'Dtrain') 
      res=rmfield(res,'Dtrain');
    end
  end
  
 case 'matlab_classify'
  tr=Dtest.supdat(supid,:)';
  [pred,res.err,res.posterior,res.logp]=...
      classify(Dtest.dat',Dtrain.dat',Dtrain.supdat(supid,:)',classifier_params.dist_type,classifier_params.prior);
  Dtest=add_D_sup(Dtest,...
                  ['PRED-' deblank(Dtest.supacc(supid,:))],...
                  ['Predicted ' deblank(Dtest.supdesc(supid,:))],...
                  pred');
  
 case 'naive_bayes'
  tr=Dtest.supdat(supid,:)';
  disp('NOT FUNCTIONAL. Use matlab_classify');
  [pred,res.score,res.sig,res.logscore]= ...
      naive_bayes_classify(Dtest.dat',Dtrain.dat',[Dtrain.supdat(supid,:); 1-Dtrain.supdat(supid,:)],...
                           classifier_params.prior,tr);
  Dtest=add_D_sup(Dtest,...
                  ['PRED-' deblank(Dtest.supacc(supid,:))],...
                  ['Predicted ' deblank(Dtest.supdesc(supid,:))],...
                  pred');
   
 case 'sp_knn'
   tr=Dtest.supdat(supid,:)';
   [cls.r,cls.a]=train(knn({distance('euclid'),'k=3'}),...
               data('train',Dtrain.dat(1: ...
                                       classifier_params.num_features,:)',2*(Dtrain.supdat(supid,:)'-0.5)));
   res=test(cls.a,data('test',Dtest.dat(1: ...
                                    classifier_params.num_features,:)',2*(tr-0.5)));
   pred=double(get_x(res)>0);
   cls.params=classifier_params;
   Dtest=add_D_sup(Dtest,...
                   ['PRED-' deblank(Dtest.supacc(supid,:))],...
                   ['Predicted ' deblank(Dtest.supdesc(supid,:))],...
                   pred');

  
 case 'gg_knn'
  class_types=unique(Dtrain.supdat(supid,:));
  if ~isempty(setdiff(class_types,[0 1]))
    [typeacc,typedesc,Dtrain,range,non_empty]= ...
        decollapse_supdat(Dtrain,supid);
    Dtrain=reorder_D_sup(Dtrain,'cols',range);
    [typeacc_ts,typedesc_ts,Dtest,range_ts,non_empty_ts]= ...
        decollapse_supdat(Dtest,1,non_empty);
    Dtest=reorder_D_sup(Dtest,'cols',range_ts);
    supid=1:length(range);
  end 
 
  tr=Dtest.supdat(supid,:)';
  if length(supid)==1
    [c,sc]=knn_classify(Dtest.dat(1:classifier_params.num_features,:)',Dtrain.dat(1:classifier_params.num_features,:)',...
                   [Dtrain.supdat(supid,:); 1-Dtrain.supdat(supid,:)],...
                   ones(2,1),classifier_params.num_neighbors,classifier_params.dist_type);
    pred=double(c==1)';
  else
    [c,cls.score,cls.knnind]=knn_classify(Dtest.dat(1:classifier_params.num_features,:)',Dtrain.dat(1:classifier_params.num_features,:)',...
                   Dtrain.supdat(supid,:),...
                   ones(length(supid),1), ...
                   classifier_params.num_neighbors,classifier_params.dist_type);
    pred=zeros(size(tr,2),length(c));
    pred(c+size(tr,2)*(0:(length(c)-1)))=1;
    pred=pred';
    [dum,ctr]=max(tr,[],2);
    nerr=0; %%%% FIXME !!!!!
  end    
  cls.params=classifier_params;
  for i=1:length(supid)
    Dtest=add_D_sup(Dtest,...
                    ['PRED-' deblank(Dtest.supacc(supid(i),:)) '(' num2str(nerr) ')'],...
                    ['Predicted ' deblank(Dtest.supdesc(supid(i),:)) ' (ERR=' num2str(nerr) ')' ],...
                    pred(:,i)');
  end
 

  
 case 'gg_pnn_w_sel_fet'
   if isfield(classifier_params,'params')
     classifier_params=add_struct(classifier_params, ...
                                  classifier_params.params);
   end
   if ischar(classifier_params.num_features)
     classifier_params.num_features=str2num(classifier_params.num_features);
   end
   tr=Dtest.supdat(supid,:)';
   cls.si=gg_sel_fet(Dtrain,supid, ...
                     classifier_params.num_features,classifier_params.fet_sel);
   Dtrain=reorder_D_rows(Dtrain,cls.si);
   Dtest=reorder_D_rows(Dtest,cls.si);
   if isfield(classifier_params,'verbose') && classifier_params.verbose
     disp('train using:');
     disp(strvcat(Dtrain.gacc));
   end
   if ~isfield(classifier_params,'class_weights')
     classifier_params.class_weights=[0.5; 0.5];
   end
   if isfield(classifier_params,'write_gct') && ...
         classifier_params.write_gct
     saveDtest=reorder_D_rows(Dtest,1: ...
                              classifier_params.num_features);
     saveDtrain=reorder_D_rows(Dtrain,1:classifier_params.num_features);
     write_mit_gct_file(['Train_' num2str(supid) '.gct'], ...
                        saveDtrain);
     write_mit_cls_file(['Train_' num2str(supid) '.cls'], saveDtrain, ...
                        supid);
     write_mit_gct_file(['Test_' num2str(supid) '.gct'],saveDtest);     
     write_mit_cls_file(['Test_' num2str(supid) '.cls'], saveDtest, ...
                        supid);
     train_test_data.Dtrain=saveDtrain;
     train_test_data.Dtest=saveDtest;
     train_test_data.supid=supid;
     train_test_data.classifier_params=classifier_params;
     save(['Data_' num2str(supid) '.mat'],'train_test_data');
   end
   [c,cls.score,cls.sig,cls.logscore]=pnn_classify(Dtest.dat(1:classifier_params.num_features,:)',Dtrain.dat(1:classifier_params.num_features,:)',...
                                                   [Dtrain.supdat(supid,:); 1-Dtrain.supdat(supid,:)],...
                                                   classifier_params.class_weights,classifier_params.sig_type,classifier_params.sig,classifier_params.dist_type,...
                                                   [Dtest.supdat(supid,:); 1-Dtest.supdat(supid,:)]);
%    disp([ supid c cls.score cls.sig cls.logscore]);
   pred=double(c==1)';
   cls.params=classifier_params;
   res=[];
   nerr=nnz(tr~=pred);
   Dtest=add_D_sup(Dtest,...
                   ['PRED-' deblank(Dtest.supacc(supid,:)) '(' num2str(nerr) ')'],...
                   ['Predicted ' deblank(Dtest.supdesc(supid,:)) ' (ERR=' num2str(nerr) ')' ],...
                   pred');
 
 case 'gg_knn_w_sel_fet'
   if isfield(classifier_params,'params')
     classifier_params=add_struct(classifier_params, ...
                                  classifier_params.params);
   end
   if ischar(classifier_params.num_features)
     classifier_params.num_features=str2num(classifier_params.num_features);
   end
   tr=Dtest.supdat(supid,:)';
   cls.si=gg_sel_fet(Dtrain,supid,classifier_params.num_features,classifier_params.fet_sel);
   Dtrain=reorder_D_rows(Dtrain,cls.si);
   Dtest=reorder_D_rows(Dtest,cls.si);
   if isfield(classifier_params,'verbose') && classifier_params.verbose
     disp('train using:');
     disp(strvcat(Dtrain.gacc));
   end
   if ~isfield(classifier_params,'class_threshold')
       classifier_params.class_threshold=classifier_params.num_neighbors/2;
   end
   w=[1; (classifier_params.class_threshold)/(classifier_params.num_neighbors-classifier_params.class_threshold)-eps]; 
   [c,cls.score,cls.knnind]=knn_classify(Dtest.dat(1:classifier_params.num_features,:)',Dtrain.dat(1:classifier_params.num_features,:)',...
                  [Dtrain.supdat(supid,:); 1-Dtrain.supdat(supid,:)],...
                  w,classifier_params.num_neighbors,classifier_params.dist_type);
   pred=double(c==1)';
   cls.params=classifier_params;
   res=[];
   nerr=nnz(tr~=pred);
   Dtest=add_D_sup(Dtest,...
                   ['PRED-' deblank(Dtest.supacc(supid,:)) '(' num2str(nerr) ')'],...
                   ['Predicted ' deblank(Dtest.supdesc(supid,:)) ' (ERR=' num2str(nerr) ')' ],...
                   pred');
  
 case 'gp_classifier'
  global my_gp
  if ~isempty(my_gp)
    % write res files
    % write cls files
    % run analysis
    % read results
    class_types=unique(Dtrain.supdat(supid,:));
    if ~isempty(setdiff(class_types,[0 1]))
      [typeacc,typedesc,Dtrain,range,non_empty]= ...
          decollapse_supdat(Dtrain,supid);
      Dtrain=reorder_D_sup(Dtrain,'cols',range);
      [typeacc_ts,typedesc_ts,Dtest,range_ts,non_empty_ts]= ...
          decollapse_supdat(Dtest,1,non_empty);
      Dtest=reorder_D_sup(Dtest,'cols',range_ts);
      supid=1:length(range);
    end 
    
    write_mit_res_file('temp_train.res',Dtrain);
    write_mit_cls_file('temp_train.cls',Dtrain,supid);
    
    write_mit_res_file('temp_test.res',Dtest);
    write_mit_cls_file('temp_test.cls',Dtest,supid);
    
    gp_params=cell2struct({'temp_train.res','temp_train.cls','temp_test.res','temp_test.cls'},...
                        {'train_filename','train_class_filename', ...
                        'test_filename','test_class_filename'},2);
    if isfield(classifier_params,'params')
      gp_params=add_struct(gp_params,classifier_params.params);
    end
    res=runAnalysis(my_gp,classifier_params.gp_classifier,gp_params);
    cls.params=classifier_params;
    [pred,tr,res]=read_mit_pred_results('temp_test.pred.odf');
    check_file_exists_and_delete('temp_train.res');
    check_file_exists_and_delete('temp_train.cls');
    check_file_exists_and_delete('temp_test.res');
    check_file_exists_and_delete('temp_test.cls');
    check_file_exists_and_delete('temp_test.pred.odf');
    nerr=nnz(tr~=pred);
    Dtest=add_D_sup(Dtest,...
                    ['PRED-' deblank(Dtest.supacc(supid,:)) '(' num2str(nerr) ')'],...
                    ['Predicted ' deblank(Dtest.supdesc(supid,:)) ' (ERR=' num2str(nerr) ')' ],...
                    pred');
  else
    error('no my_gp');
  end
  
  
%  case 'gp_knn'
%   global my_gp
%   if ~isempty(my_gp)
%     % write res files
%     % write cls files
%     % run analysis
%     % read results
%     write_mit_res_file('temp_train.res',Dtrain);
%     write_mit_cls_file('temp_train.cls',Dtrain,supid);
    
%     write_mit_res_file('temp_test.res',Dtest);
%     write_mit_cls_file('temp_test.cls',Dtest,supid);
    
%     gp_params=cell2struct({'temp_train.res','temp_train.cls','temp_test.res','temp_test.cls'},...
%                         {'train_filename','train_class_filename', ...
%                         'test_filename','test_class_filename'},2);
%     if isfield(classifier_params,'params')
%       gp_params=add_struct(gp_params,classifier_params.params);
%     end
%     res=runAnalysis(my_gp,'KNN',gp_params);
%     cls.params=classifier_params;
%     [pred,tr,res]=read_mit_pred_results('temp_test.pred.odf');
%     check_file_exists_and_delete('temp_train.res');
%     check_file_exists_and_delete('temp_train.cls');
%     check_file_exists_and_delete('temp_test.res');
%     check_file_exists_and_delete('temp_test.cls');
%     check_file_exists_and_delete('temp_test.pred.odf');
%   else
%     error('no my_gp');
%  end
 
  
 otherwise 
  cls=train(Dtrain,supid,classifier_params);
  [pred,tr,res]=predict(Dtest,classifier,supid);
end



