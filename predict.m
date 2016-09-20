function [pred,tr,res]=predict(Dtest,classifier,supid)

if nnz(isnan(D.supdat(supid,:)))
  error('does not support NaNs');
end
ns=size(D.supdat,2);

if ischar(classifier_params)
  classifier_params.method=classifier_params;
end

switch classifier.classifier_params.method
 case 'cv_select_features'
  Dtest_f=reorder_D(Dtest,'rows',classifier.use_features);
  [pred,tr,res]=predict(Dtest_f,classifier.use_classifier,supid);
 case 'gp_knn'
  error('must run train_and_predict');
 otherwise
  error('no such prediction method');
end
