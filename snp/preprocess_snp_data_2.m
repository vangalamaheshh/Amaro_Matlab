function [CL,P,batches,snps,all_mn,notref]=preprocess_snp_data_2(M,perform_correct_batch_effect,divide_by_normals)

if ~exist('perform_correct_batch_effect','var');
  perform_correct_batch_effect=1;
end

if ~exist('divide_by_normals','var');
  divide_by_normals=1;
end

if iscell(M)
  clear M;
  global M;
  use_global=1;
%  global gM;
%  ML1=gM{M{1}};
else
  use_global=0;
end

if median(M.dat(:))>0.5
  M.dat(M.dat<1)=1;
  M.dat=log2(M.dat)-1;
end

% multiply by 2 (add 1 in log2) for males on X chromosome
% FIXME: don't add 1 but the difference between the median of
% median on X chormosomes between the males and females.
Xpos=find(M.chrn==23);
M.dat(Xpos,:)=M.dat(Xpos,:)+repmat(M.supdat(find_supid(M,'GENDER'),:)==1,length(Xpos),1); % Gender:
                                                                                                  % male=1, female=2 

if (0) % scale samples 
  x=median(M.dat,1);
  mx=median(x);
  M.dat=M.dat-repmat(x-mx,size(M.dat,1),1);
end

if perform_correct_batch_effect
  disp('correcting batch effect');
  if use_global
    [M,P,batches,snps]=correct_batch_effect({},5,0.05,0.001);
  else
    [M,P,batches,snps]=correct_batch_effect(ML2,5,0.05,0.001);
  end
else
  P=[];
  batches=[];
  snps=[];
end

% divide by normals
ref=find(M.supdat(find_supid(M,'N'),:));
notref=setdiff(1:size(M.dat,2),ref);

REF=reorder_D_cols(M,ref);
CL=reorder_D_cols(M,notref);

if divide_by_normals
  all_mn=mean(REF.dat,2);
  uc=unique(CL.supdat(find_supid(CL,'CORE'),:));
  for i=1:length(uc)
    %%% FIXME: add inner loop on tissue type
    cl_incore=find(CL.supdat(find_supid(CL,'CORE'),:)==uc(i));
    ref_incore=find(REF.supdat(find_supid(REF,'CORE'),:)==uc(i));
    if isempty(ref_incore)
      CL.dat(:,cl_incore)=CL.dat(:,cl_incore)-repmat(all_mn,1,length(cl_incore));
      disp('Could not find reference for samples:');
      disp(catn(CL.sdesc(cl_incore)));
    else
      disp([ 'Core ' num2str(uc(i)) ') Dividing ' num2str(length(cl_incore)) ' by ' num2str(length(ref_incore)) ' reference samples']);
      CL.dat(:,cl_incore)=CL.dat(:,cl_incore)-repmat(mean(REF.dat(:,ref_incore),2),1,length(cl_incore));
    end
  end
else
  disp('Not dividing by normals');
  all_mn=[];
end
% previously was:
% refdat(Xpos,:)=refdat(Xpos,:)+repmat(M.supdat(find_supid(M,'GENDER'),ref),length(Xpos),1); 
% 
% CL.dat=CL.dat-repmat(mn,1,size(CL.dat,2)); % subtract since in log space
  
CL.dat=2.^(CL.dat+1); % back to copy number





