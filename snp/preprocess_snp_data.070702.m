function [CL,P,batches,snps,all_mn,notref]=preprocess_snp_data(M,params); 
%
%perform_correct_batch_effect,divide_by_normals,use_median_scale,use_n_close_normals,use_matched_normal,use_all_cores_for_ref,use_median_of_ref,use_tangent_normalization)
%


if ~exist('perform_correct_batch_effect','var')
  perform_correct_batch_effect=1;
end

if ~exist('divide_by_normals','var')
  divide_by_normals=1;
end

if ~exist('use_tangent_normalization','var')
  use_tangent_normalization=0;
end

if ~exist('use_median_scale','var')
  use_median_scale=0;
end

if ~exist('use_n_close_normals','var')
  use_n_close_normals=0;
end

if ~exist('use_matched_normal','var')
  use_matched_normal=0;
end

if ~exist('use_all_cores_for_ref','var')
  use_all_cores_for_ref=0;
end

if ~exist('use_median_of_ref','var')
  use_median_of_ref=0;
end

if ~exist('use_tangent_normalization','var')
  use_tangent_normlization=0;
end

if iscell(M)
  clear M;
  global M;
  use_global=1;
  ML1=M;
else
  ML1=M;
  clear M;
  use_global=0;
end

% assuming M contains signal (M>>1)
ML1.dat(ML1.dat<1)=1; % to avoid small numbers 
ML1.dat=log2(ML1.dat);

% multiply by 2 (add 1 in log2) for males on X chromosome
% FIXME: don't add 1 but the difference between the median of
% median on X chormosomes between the males and females.
disp('adding 1 to the X chrom. of males');
Xpos=find(ML1.chrn==23);
ML1.dat(Xpos,:)=ML1.dat(Xpos,:)+repmat(ML1.supdat(find_supid(ML1,'GENDER'),:)==1,length(Xpos),1); % Gender:
                                                                                                  % male=1, female=2 
% Subtract median before batch correction
% rescaling the scaling
% overrides brightness correction
if params.median_scale_before_batch_correct
  disp('subtracting median of values');
  x=median(ML1.dat,1);
%  mx=median(x);
  ML2=ML1;
  ML2.dat=ML2.dat-repmat(x,size(ML2.dat,1),1);
else
  ML2=ML1;
end
clear ML1;

if params.perform_batch_correction
  disp('correcting batch effect');
  if use_global
    M=ML2;
    clear ML2;
    [ML3,P,batches,snps]=correct_batch_effect({},5,0.05,0.001);    
  else
    [ML3,P,batches,snps]=correct_batch_effect(ML2,5,0.05,0.001);
  end
else
  ML3=ML2;
  P=[];
  batches=[];
  snps=[];
end

if (0)
  save('ML3.mat','ML3');
end

clear ML2;

% divide by normals
if use_n_close_normals==0
  ref=find(ML3.supdat(find_supid(ML3,'N'),:));
  notref=setdiff(1:size(ML3.dat,2),ref);

  REF=reorder_D_cols(ML3,ref);
  CL=reorder_D_cols(ML3,notref);
else
  if isfield(ML3.sis,'cellxeno')
    ref=intersect(find(ML3.supdat(find_supid(ML3,'N'),:) | ML3.supdat(find_supid(ML3,'CTRL'),:)),grep('no', ...
                                                      {ML3.sis.cellxeno},1));
  else
    ref=find(ML3.supdat(find_supid(ML3,'N'),:) | ML3.supdat(find_supid(ML3,'CTRL'),:));
  end
  notref=setdiff(1:size(ML3.dat,2),ref);
  REF=reorder_D_cols(ML3,ref);
  CL=ML3;
end

if divide_by_normals
  all_mn=mean(REF.dat,2);
  uc=unique(CL.supdat(find_supid(CL,'CORE'),:));
  if use_all_cores_for_ref
    uc=1;
  end
  for i=1:length(uc)
    %%% FIXME: add inner loop on tissue type
    if use_all_cores_for_ref
      cl_incore=1:size(CL.dat,2);
      ref_incore=1:size(REF.dat,2);
    else
      cl_incore=find(CL.supdat(find_supid(CL,'CORE'),:)==uc(i));
      ref_incore=find(REF.supdat(find_supid(REF,'CORE'),:)==uc(i));
    end
    if isempty(ref_incore)
      CL.dat(:,cl_incore)=CL.dat(:,cl_incore)-repmat(all_mn,1,length(cl_incore));
      disp('Could not find reference for samples:');
      disp(catn(CL.sdesc(cl_incore)));
    else
      if use_n_close_normals>0
        d=dist(REF.dat(CL.chrn~=23,ref_incore)',CL.dat(CL.chrn~=23,cl_incore)','euclidean_after_med_subtract');
        d(d==0)=Inf; % dont use itself to normalize
        [sd,si]=sort(d,1);
        for i=1:length(cl_incore)
          if use_matched_normal 
            % find matched normal
            use_normals=strmatch(CL.sis(cl_incore(i)).paired,REG.sdesc,'exact'); % fix field name matched...
          end
          if ~use_matched_normal || isempty(use_normals)
            use_normals=ref_incore(si(1:min(size(si,1),use_n_close_normals),i));
            %          d1=dist(REF.dat(CL.chrn~=23,use_normals)',CL.dat(CL.chrn~=23,cl_incore(i))', ...
            %                  'euclidean_after_med_subtract');
          end
          if use_median_of_ref
            CL.dat(:,cl_incore(i))=CL.dat(:,cl_incore(i))-median(REF.dat(:,use_normals),2);
          else
            CL.dat(:,cl_incore(i))=CL.dat(:,cl_incore(i))-mean(REF.dat(:,use_normals),2);
          end
          CL.used_normals{i}=REF.sdesc(use_normals);
          CL.temp{i,1}=REF.sdesc(ref_incore(si(1:min(size(si,1),Inf),i)));
          CL.temp{i,2}=sd(:,i);          
        end
      else
        if use_tangent_normalization
          disp('using tangent normalization');
          xx=tangent_normalization(CL.dat(:,cl_incore),REF.dat(:,ref_incore));
          CL.dat(:,cl_incore)=xx;
        else
          disp([ 'Core ' num2str(uc(i)) ') Dividing ' num2str(length(cl_incore)) ' by ' num2str(length(ref_incore)) ' reference samples']);
          if use_median_of_ref
            CL.dat(:,cl_incore)=CL.dat(:,cl_incore)-repmat(median(REF.dat(:,ref_incore),2),1,length(cl_incore));
          else
            CL.dat(:,cl_incore)=CL.dat(:,cl_incore)-repmat(mean(REF.dat(:,ref_incore),2),1,length(cl_incore));
          end
        end
      end
    end
  end
else
  disp('Not dividing by normals');
  all_mn=[];
end
% previously was:
% refdat(Xpos,:)=refdat(Xpos,:)+repmat(ML3.supdat(find_supid(ML3,'GENDER'),ref),length(Xpos),1); 
% 
% CL.dat=CL.dat-repmat(mn,1,size(CL.dat,2)); % subtract since in log space
  
if params.median_scale_after_norm
  disp('subtracting median of values (after dividing by normals) excluding the X');
  CL.dat=CL.dat-repmat(median(CL.dat(CL.chrn~=23,:),1),size(CL.dat,1),1);
end

CL.dat=2.^(CL.dat+1); % back to copy number




