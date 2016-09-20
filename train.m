function cls=train(D,supid,classifier_params)

if ischar(classifier_params)
  classifier_params.method=classifier_params;
end

cls.params=classifier_params;
switch classifier_params.method
 case 'cv_select_features'
  cls0=find(D.supdat(supid,:)==0);
  cls1=find(D.supdat(supid,:)==1);
  [P1sgte,P1slte,P2s,Pf,rs,Pf2,S,gp,fwer,fpr]=marker_selection(D.dat,cls0,cls1,...
                                                    classifier_params.marker_selection.test_type,...
                                                    classifier_params.marker_selection.nperm)
  [P2s_sort,feat_ord]=sort(P2s(:,1));
  fres=classifier_params.f_start:classifier_params.f_step: ...
       min(classifier_params.f_end,length(feat_ord));
  
  for f=1:size(fres,2)
    Df=reorder_D(D,'rows',feat_ord(1:f));
    [f_confmat,f_nerr,f_cvres]=crossvalidate(Df,supid, ...
                                             classifier_params.cv_params,classifier_params.classifier_params);
    fres(2,f)=f_nerr;
  end
  f_choose_i=min(fres(2,:));
  cls.feat_ord=feat_ord;
  cls.fres=fres;
  cls.f_n_choosen=fres(1,f_choose_i);
  cls.use_features=feat_ord(1:cls.f_n_choosen);
%  cls.use_classifier
 case 'gp_cv_select_features'
  error('use train_and_predict');
 case 'gp_knn'
  error('use train_and_predict');
 case 'gp_weighted_voting'
  error('use train_and_predict');
 otherwise 
  error('none!');
end  
  
