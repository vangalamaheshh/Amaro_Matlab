function [pv,dg,class,bonf_pv] = differential_region_analysis_permutations(D,class_idx,nperm)
  
  control_idx = setdiff(1:size(D.dat,2),class_idx);
  
  if isempty(control_idx)
    error('Must supply control samples!');
  end
  
  n1 = length(control_idx);
  n2 = length(class_idx);
  n = n1+n2;
  
  m0=nanmean(D.dat(:,control_idx),2);
  m1=nanmean(D.dat(:,class_idx),2);
  dg = m1-m0;
  
  pv = zeros(1,size(D.dat,1));
  
  perm_dgs = zeros(length(dg),nperm);
  
  disp(['Permutation: 0 of ' num2str(nperm)])
  for i=1:nperm
    if mod(i,1000) == 0
      disp(['Permutation: ' num2str(i) ' of ' num2str(nperm)])
    end
    cur_perm = randperm(n);
    cur_m0 = nanmean(D.dat(:,cur_perm(1:n1)),2);
    cur_m1 = nanmean(D.dat(:,cur_perm((n1+1):end)),2);
    perm_dgs(:,i) = cur_m1-cur_m0;
  end
  
  for j=1:length(dg)
    qq = find(abs(perm_dgs(j,:)));
    if ~isempty(qq)
      pv(j) = length(find(abs(perm_dgs(j,:)) > abs(dg(j))))/nperm;
    else
      pv(j) = 1;
    end
    
    if pv(j) == 0
      pv(j) = 1/nperm;
    end
    bonf_pv(j) = min(1,pv(j)*size(D.dat,1));
    if sign(dg(j))~=0
      class(j) = (3-sign(dg(j)))/2;
    else
      class(j) = 3;
    end
    
  end
  
    
