function [train_idx,test_idx]=get_train_test_split(D,supid,params)

switch params.method
 case 'test_n_class'
  v=params.n_class_vec;
  idx=[];
  for i=1:length(v)
    if v(i)>0
      class_idx=find(D.supdat(supid,:)==i);
      if ~isempty(class_idx)
        r=randperm(length(class_idx));
        idx=[ idx class_idx(r(1:v(i))) ];
      end
    end
  end
  test_idx=idx;
  train_idx=setdiff(1:size(D.dat,2),test_idx);
 case 'test_fraction'
  v=histc(D.supdat(supid,:),1:max(D.supdat(supid,:)));
  vtest=floor(v*params.test_fraction);
  if isfield(params,'min_test')
    below_min=find(vtest<params.min_test);
    if ~isempty(below_min)
      for i=1:length(below_min)
        if v(below_min(i))<=params.min_test
          verbose(['Not enough samples in type #' num2str(i)]);
        else
          vtest(below_min(i))=params.min_test;
        end
      end
    end
  end
  [train_idx,test_idx]=get_train_test_split(D,supid,struct('method','test_n_class','n_class_vec',vtest));
 otherwise
  error('no such method');
end

