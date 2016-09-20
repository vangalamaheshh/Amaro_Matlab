function [confmat,nerr,cvres,D,params]=crossvalidate(D,supid,cv_params)

if nnz(isnan(D.supdat(supid,:)))
  error('does not support NaNs');
end
ns=size(D.supdat,2);

params=[];
% assume LOOCV
if ischar(cv_params)
  cv_params.method=cv_params;
end
class_types=unique(D.supdat(supid,:));
for i=1:length(class_types)
  class_enums{i,1}=class_types(i);
  class_enums{i,2}=i;
end

ntypes=length(class_types);
switch cv_params.method
 case 'loocv'
  confmat=zeros(ntypes,ntypes);
  if isfield(cv_params,'lsfdir')
    l=lsf(cv_params.lsfdir,[],0); % 0= do not use compiler
    l=maxjobs(l,cv_params.lsf_maxjobs);
    h=zeros(ns,1);
    for i=1:ns
      idx=setdiff(1:ns,i);
      Dtrain=reorder_D(D,'cols',idx);
      Dtrain.sel_fet_id=i; 
      Dtest=reorder_D(D,'cols',i);
      % HACK!!! sending all of D and picking the result for i
      if isfield(cv_params,'params')
        cv_params.classifier=add_struct(cv_params.classifier, ...
                                        cv_params.params);
      end
      if strcmp(cv_params.classifier.method,'gp_classifier')
        [l,h(i)]=bsub(l,{'pred','tr','res','classifier'}, ...
                      'train_and_predict',{Dtrain,D,supid, ...
                            cv_params.classifier});
      else
        [l,h(i)]=bsub(l,{'pred','tr','res','classifier'}, ...
                      'train_and_predict',{Dtrain,Dtest,supid, ...
                            cv_params.classifier});
      end
    end
    [l,lres]=wait(l);
    for i=1:ns
      if strcmp(cv_params.classifier.method,'gp_classifier')
        cvres(i).pred=lres{h(i)}.pred(i);
        cvres(i).true=lres{h(i)}.tr(i);
      else
        cvres(i).pred=lres{h(i)}.pred;
        cvres(i).true=lres{h(i)}.tr;        
      end
      cvres(i).predres=lres{h(i)}.res;
      cvres(i).classifier=lres{h(i)}.classifier;
      %    disp([i cvres(i).pred cvres(i).true]);
      confmat=confmat+crosstab_full(convert_enum(cvres(i).pred, ...
                                                 class_enums), ...
                                    convert_enum(cvres(i).true, ...
                                                 class_enums),1:ntypes);
      if isfield(cv_params,'save_every_n_iter') && ...
            mod(i,cv_params.save_every_n_iter)==0
        save(['current_state_' num2str(i) '.mat'],'cvres','confmat');
      end
    end   
  else
    for i=1:ns
      idx=setdiff(1:ns,i);
      Dtrain=reorder_D(D,'cols',idx);
      Dtrain.sel_fet_id=i;
      Dtest=reorder_D(D,'cols',i);
      % HACK!!! sending all of D and picking the result for i
      if isfield(cv_params,'params')
        cv_params.classifier=add_struct(cv_params.classifier, ...
                                        cv_params.params);
      end
      if strcmp(cv_params.classifier.method,'gp_classifier')
        [pred,tr,res,classifier]=train_and_predict(Dtrain,D,supid,cv_params.classifier);
        cvres(i).pred=pred(i);
        cvres(i).true=tr(i);
      else
        [pred,tr,res,classifier]=train_and_predict(Dtrain,Dtest,supid,cv_params.classifier);
        cvres(i).pred=pred;
        cvres(i).true=tr;
      end      
      cvres(i).predres=res;
      cvres(i).classifier=classifier;
      %    disp([i cvres(i).pred cvres(i).true]);
      
      fprintf(1,'.');
      confmat=confmat+crosstab_full(convert_enum(cvres(i).pred,class_enums), ...
                                    convert_enum(cvres(i).true,class_enums),1:ntypes);
      
      if isfield(cv_params,'save_every_n_iter') && ...
            mod(i,cv_params.save_every_n_iter)==0
        save(['current_state_' num2str(i) '.mat'],'cvres','confmat');
      end
    end
  end
  nerr=sum(confmat(:))-sum(diag(confmat));
  D=add_D_sup(D,['CV(' num2str(nerr) ')'],['Cross-validation (ERR=' ...
                      num2str(nerr) ')'],cat(2,cvres(:).pred));

 case 'train_test_split'
  n_runs=cv_params.n_runs;
  confmat=zeros(ntypes,ntypes);
  if isfield(cv_params,'lsfdir')
    l=lsf(cv_params.lsfdir,[],0); % 0= do not use compiler
    l=maxjobs(l,cv_params.lsf_maxjobs);
    h=zeros(n_runs,1);
    for i=1:n_runs
      [train_idx,test_idx]=get_train_test_split(D,supid,cv_params.split);
      Dtrain=reorder_D(D,'cols',train_idx);
      Dtest=reorder_D(D,'cols',test_idx);
      % HACK!!! for GP sending all of D and picking the result for test_idx
      if isfield(cv_params,'params')
        cv_params.classifier=add_struct(cv_params.classifier, ...
                                        cv_params.params);
      end
      if strcmp(cv_params.classifier.method,'gp_classifier')
        [l,h(i)]=bsub(l,{'pred','tr','res','classifier'}, ...
                      'train_and_predict',{Dtrain,D,supid, ...
                            cv_params.classifier});
      else
        [l,h(i)]=bsub(l,{'pred','tr','res','classifier'}, ...
                      'train_and_predict',{Dtrain,Dtest,supid, ...
                            cv_params.classifier});
      end
    end
    [l,lres]=wait(l);
    for i=1:n_runs
      if strcmp(cv_params.classifier.method,'gp_classifier')
        cvres(i).pred=lres{h(i)}.pred(test_idx);
        cvres(i).true=lres{h(i)}.tr(test_idx);
      else
        cvres(i).pred=lres{h(i)}.pred;
        cvres(i).true=lres{h(i)}.tr;        
      end
      cvres(i).predres=lres{h(i)}.res;
      cvres(i).classifier=lres{h(i)}.classifier;
      %    disp([i cvres(i).pred cvres(i).true]);
      [tmp,pred_ind]=max(cvres(i).pred,[],2);
      [tmp,tr_ind]=max(cvres(i).true,[],2);
      cv_confmat=crosstab_full(convert_enum(pred_ind,class_enums), ...
                               convert_enum(tr_ind,class_enums),1:ntypes);
      cvres(i).confmat=cv_confmat;
      cvres(i).nerr=sum(cv_confmat(:))-sum(diag(cv_confmat));
      confmat=confmat+cv_confmat;
      if isfield(cv_params,'save_every_n_iter') && ...
            mod(i,cv_params.save_every_n_iter)==0
        save(['current_state_' num2str(i) '.mat'],'cvres','confmat');
      end
    end   
  else
    for i=1:n_runs
      [train_idx,test_idx]=get_train_test_split(D,supid,cv_params.split);
      Dtrain=reorder_D(D,'cols',train_idx);
      Dtest=reorder_D(D,'cols',test_idx);
      cvres(i).test_idx=test_idx;
      cvres(i).train_idx=train_idx;
      % HACK!!! for GP sending all of D and picking the result for test_idx
      if isfield(cv_params,'params')
        cv_params.classifier=add_struct(cv_params.classifier, ...
                                        cv_params.params);
      end
      if strcmp(cv_params.classifier.method,'gp_classifier')
        [pred,tr,res,classifier]=train_and_predict(Dtrain,D,supid,cv_params.classifier);
        cvres(i).pred=pred(i);
        cvres(i).true=tr(i);
      else
        [pred,tr,res,classifier]=train_and_predict(Dtrain,Dtest,supid,cv_params.classifier);
        cvres(i).pred=pred;
        cvres(i).true=tr;
      end      
      cvres(i).predres=res;
      if ~isfield(cv_params,'dont_save_classifier') || ~cv_params.dont_save_classifier
        cvres(i).classifier=classifier;
      end
      
      fprintf(1,'.');
      [tmp,pred_ind]=max(cvres(i).pred,[],2);
      [tmp,tr_ind]=max(cvres(i).true,[],2);
      cv_confmat=crosstab_full(convert_enum(pred_ind,class_enums), ...
                               convert_enum(tr_ind,class_enums),1:ntypes);
      cvres(i).confmat=cv_confmat;
      cvres(i).nerr=sum(cv_confmat(:))-sum(diag(cv_confmat));
      confmat=confmat+cv_confmat;
      if isfield(cv_params,'save_every_n_iter') && ...
            mod(i,cv_params.save_every_n_iter)==0
        save(['current_state_' num2str(i) '.mat'],'cvres','confmat');
      end
    end   
  end
  nerr=sum(confmat(:))-sum(diag(confmat));
  D=add_D_sup(D,['CV(' num2str(nerr) ')'],['Cross-validation (ERR=' ...
                      num2str(nerr) ')'],cat(2,cvres(:).pred));

% case 'gp_xv' %% FIX ME: make a genertic gp_xv
 case 'gp_weighted_voting_xv'
  confmat=zeros(ntypes,ntypes);
  gp_params=struct('data_filename','temp.res','class_filename','temp.cls');
  if ~isfield(cv_params,'params') | ...
        ( isfield(cv_params,'params') && ...
          ( ~isfield(cv_params.params,'data_filename') | ...
            ~isfield(cv_params.params,'class_filename'))) 
    write_mit_res_file('temp.res',D);
    write_mit_cls_file('temp.cls',D,supid);
    if isfield(cv_params,'presend_files') && ...
          cv_params.presend_files
      gp_params=gp_presend_files(gp_params);
      params=gp_params;
    end
  end
  if isfield(cv_params,'params')
    gp_params=add_struct(gp_params,cv_params.params);
  end
  global my_gp;
  if ~isempty(my_gp)
    res=runAnalysis(my_gp,'WeightedVotingXValidation',gp_params);
    
    [pred,tr,res]=read_mit_pred_results('temp.pred.odf');
    confmat=crosstab_full(pred+1,tr+1,1:2);
    nerr=sum(confmat(:))-sum(diag(confmat));
    for i=1:length(pred)
      cvres(i).pred=pred(i);
      cvres(i).true=tr(i);
    end;
    D=add_D_sup(D,['CV(' num2str(nerr) ')'],['Cross-validation (ERR=' ...
                        num2str(nerr) ')'],cat(2,cvres(:).pred));
  else
    error('no my_gp');
  end

 case 'gp_knn_xv'
  confmat=zeros(ntypes,ntypes);
  gp_params=struct('data_filename','temp.res','class_filename','temp.cls');
  if ~isfield(cv_params,'params') | ...
        ( isfield(cv_params,'params') && ...
          ( ~isfield(cv_params.params,'data_filename') | ...
            ~isfield(cv_params.params,'class_filename'))) 
    write_mit_res_file('temp.res',D);
    write_mit_cls_file('temp.cls',D,supid);
    if isfield(cv_params,'presend_files') && ...
          cv_params.presend_files
      gp_params=gp_presend_files(gp_params);
      params=gp_params;
    end
  end
  if isfield(cv_params,'params')
    gp_params=add_struct(gp_params,cv_params.params);
  end
    
  global my_gp;
  if ~isempty(my_gp)
    res=runAnalysis(my_gp,'KNNXValidation',gp_params);
%    res=runAnalysis(my_gp,'KNNXValidation2',gp_params);
    [pred,tr,res]=read_mit_pred_results('temp.pred.odf');
    confmat=crosstab_full(pred+1,tr+1,1:2);
    nerr=sum(confmat(:))-sum(diag(confmat));
    for i=1:length(pred)
      cvres(i).pred=pred(i);
      cvres(i).true=tr(i);
    end;
    D=add_D_sup(D,['CV(' num2str(nerr) ')'],['Cross-validation (ERR=' ...
                        num2str(nerr) ')'],cat(2,cvres(:).pred));
  else
    error('no gp');
  end
  
end

