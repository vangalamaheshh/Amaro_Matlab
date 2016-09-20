function [C2,res]=SNP60_pipeline(H,params)

if ischar(params)
  tmp=struct('method',params);
  params=tmp;
end

switch params.method
 case 'v2'
  if ischar(params.birdseed_clusters)
    B=read_birdseed_clusters(params.birdseed_clusters);
  else
    B=params.birdseed_clusters;
  end
  res.birdseed_clusters=B;
  b=strvcat(B.marker);
  disp([ 'Using ' num2str(size(b,1)) ' birdseed clusters']);
  % unique(b(:,1)) % S
  dash2=find(b(:,end)=='2');
  disp([ 'Found ' num2str(length(dash2)) ' with -2 at the end']);
  b_num=str2int_matrix(b(dash2,1:(end-2))); % remove -2 at end
  H_snp_idx=find(H.gsupdat(1,:));
  tmp=H.gsupdat(3,H_snp_idx);
  [St,s1,s2]=match_num_sets(b_num,tmp);
  disp([ num2str(length(s2)) ' match H which has ' num2str(length(H_snp_idx)) ' SNPs']);
%   range(b_num(s1)-tmp(s2)') % 0
  C=reorder_D_rows(H,H_snp_idx(s2));
  if isfield(C,'orig')
    C=rmfield(C,{'orig'});
  end
  B1=reorder_D_rows(B,dash2(s1));
  %same(C.marker,b(s1,1:(end-2))) % 1
  
  x0=B1.dat(:,[13 2]);
  x1=B1.dat(:,[7 8]);
  x2=B1.dat(:,[1 14]);
  
  keyboard
  tmp=[x0(:,1) x1(:,1) x2(:,1)];
  [stmp,sti]=sort(tmp(:,2));
  stmp=tmp(sti,:);
  d=stmp(:,2:3)-repmat(stmp(:,1),1,2);
  r=d(:,2)./(d(:,1)+eps);
  
  disp(['replacing ' num2str(nnz(isnan(x0))) ' NaNs in the background terms with 0']);
  x0(isnan(x0))=0;
  
  % deltas 
  delta01=x1-x0;
  delta12=x2-x1;
  delta02=x2-x0;
  
  delta02_div_2=nanmean(cat(3,delta01,delta12),3); %delta02/2;
  delta02_div_2_avg=repmat(nanmean(delta02_div_2,2),1,2);

%  min_val=min(C.adat,[],
  
  bad_probes=find(any(delta02_div_2<=params.delta_cutoff,2));
  good_probes=setdiff(1:size(C.dat,1),bad_probes);

  C1=reorder_D_rows(C,good_probes);
  C1.adat=C1.adat-repmat(permute(x0(good_probes,:),[1 3 2]),[1 size(C1.adat,2) 1]);
  C1.adat=C1.adat./repmat(permute(delta02_div_2(good_probes,:),[1 3 2]),[1 size(C1.adat,2) 1]);
  disp(['Applying allele-specific floor value of ' num2str(params.floor_val)]);
  C1.adat(C1.adat<params.floor_val)=params.floor_val;
  C1.dat=sum(C1.adat,3);


  C2=C1;
  C2.dat=log2(C2.dat)-1;
  
  res.bad_probes=bad_probes;
  res.delta=delta02_div_2;
  
 otherwise
  error('not implemented');
end
