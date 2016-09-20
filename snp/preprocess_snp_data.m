function [CL,P,batches,snps,all_center,notref]=preprocess_snp_data(M,params); 
%
%perform_correct_batch_effect,divide_by_normals,use_median_scale,use_n_close_normals,use_matched_normal,use_all_cores_for_ref,use_median_of_ref,use_tangent_normalization)
%

if isfield(params,'input_is_raw')
  input_is_raw=params.input_is_raw;
else
  input_is_raw=1;
end

if isfield(params,'perform_batch_correction')
  perform_correct_batch_effect = params.perform_batch_correction;
else
  perform_correct_batch_effect=0;
end

if perform_correct_batch_effect
  if isfield(params,'batch_effect_correction')  
    batch_effect_params = params.batch_effect_correction;
  else
    error(['Must supply batch effect parameters if running batch effect ' ...
           'correction!']);
  end
end

if isfield(params,'divide_by_normals')
  divide_by_normals = params.divide_by_normals;
else 
  divide_by_normals=1;
end

if isfield(params,'collapse')
  norm_collapse_method = params.collapse;
else
  norm_collapse_method = 'median';
end

if isfield(params,'selection')
  norm_select_method = params.selection;
else
  norm_select_method = 'closest_n';
end

if isfield(params,'use_paired')
  use_paired = params.use_paired;
else
  use_paired = 0;
end

if isfield(params,'n_closest_n')
  n_closest_n = params.n_closest_n;
else
  n_closest_n = 0;
end

if isequal(norm_select_method,'closest_n') && n_closest_n < 1
  error('n_closest_n should be greater or equal to 1');
end

if isfield(params,'normalize_tumors')
  normalize_tumors=params.normalize_tumors;
else
  normalize_tumors=1;
end

if isfield(params,'normalize_normals')
  normalize_normals=params.normalize_normals;
else
  normalize_normals=0;
end
  
if isfield(params,'use_all_cores_for_ref')
  use_all_cores_for_ref = params.use_all_cores_for_ref;
else
  use_all_cores_for_ref = 1;
end

if isfield(params,'median_scale_before_batch_correct')
  median_scale_before_bc = params.median_scale_before_batch_correct;
else
  median_scale_before_bc = 1;
end

if isfield(params,'median_scale_after_norm')
  median_scale_after_norm = params.median_scale_after_norm;
else
  median_scale_after_norm = 0;
end

if isfield(params,'ignore_flagged_normals')
  ignore_flagged_normals = params.ignore_flagged_normals;
else
  ignore_flagged_normals = 1;
end

if isfield(params,'conserve_memory') && params.conserve_memory
  conserve_memory=1;
else
  conserve_memory=0;
end

% check if to use a global M
if iscell(M)
  if isempty(M)
    clear M;
    global M;
    use_global=1;
    ML1=M;
  elseif conserve_memory
    disp('Using conserve memory option');
  else
    error('Cannot process M');
  end
else
  ML1=M;
  clear M;
  use_global=0;
end

if conserve_memory
  % assuming M contains signal (M>>1)
  for pli=1:length(M{2})
    disp(M{2}{pli});
    ML1=load_D([ M{1} M{2}{pli} M{3} '.mat']);
    % assuming M contains signal (M>>1)
    if input_is_raw
      disp('Assuming input is raw signal.. replacing values <1 with 1');
      ML1.dat(ML1.dat<params.floor_val)=params.floor_val; % to avoid small numbers 
    end
    disp('Taking log2 of data');
    ML1.dat=log2(ML1.dat)-1; % since we use +1 at the end
    
    % multiply by 2 (add 1 in log2) for males on X chromosome
    % FIXME: don't add 1 but the difference between the median of
    % median on X chormosomes between the males and females.
    if input_is_raw
      disp('adding 1 to the X chrom. of males');
      Xpos=find(ML1.chrn==23);
      ML1.dat(Xpos,:)=ML1.dat(Xpos,:)+repmat(ML1.supdat(find_supid(ML1,'GENDER'),:)==1,length(Xpos),1); % Gender:
    end
    % male=1, female=2 
    % Subtract median before batch correction
    % rescaling the scaling
    % overrides brightness correction
    if perform_correct_batch_effect && median_scale_before_bc
      disp('subtracting median of values');
      x=median(ML1.dat,1);
      %  mx=median(x);
      ML2=ML1;
      ML2.dat=ML2.dat-repmat(x,size(ML2.dat,1),1);
    else
      ML2=ML1;
    end
    clear ML1;
    save_D([ M{1} M{2}{pli} '.ML2.mat'],ML2);
  end
else
  % assuming M contains signal (M>>1)
  if input_is_raw
    disp('Assuming input is raw signal.. replacing values <1 with 1');
    ML1.dat(ML1.dat<params.floor_val)=params.floor_val; % to avoid small numbers 
  end
  disp('Taking log2 of data');
  ML1.dat=log2(ML1.dat)-1; % since we use +1 at the end
  
  % multiply by 2 (add 1 in log2) for males on X chromosome
  % FIXME: don't add 1 but the difference between the median of
  % median on X chormosomes between the males and females.
  if input_is_raw
    disp('adding 1 to the X chrom. of males');
    Xpos=find(ML1.chrn==23);
    ML1.dat(Xpos,:)=ML1.dat(Xpos,:)+repmat(ML1.supdat(find_supid(ML1,'GENDER'),:)==1,length(Xpos),1); % Gender:
  end
  % male=1, female=2 
  % Subtract median before batch correction
  % rescaling the scaling
  % overrides brightness correction
  if perform_correct_batch_effect && median_scale_before_bc
    disp('subtracting median of values');
    x=median(ML1.dat,1);
    %  mx=median(x);
    ML2=ML1;
    ML2.dat=ML2.dat-repmat(x,size(ML2.dat,1),1);
  else
    ML2=ML1;
  end
  clear ML1;
end
    

if perform_correct_batch_effect
  disp('correcting batch effect');
  if use_global
    M=ML2;
    clear ML2;
    [ML3,P,batches,snps]=correct_batch_effect({},batch_effect_params);    
  elseif conserve_memory
    M{3}='.ML2';
    batch_effect_params.conserve_memory=1;
    [ML3,P,batches,snps]=correct_batch_effect(M,batch_effect_params);       
  else
    [ML3,P,batches,snps]=correct_batch_effect(ML2,batch_effect_params);
  end
else
  ML3=ML2;
  P=[];
  batches=[];
  snps=[];
end

save_D('ML3.mat',ML3,'-v7.3');
keyboard

clear ML2;

% divide by normals
if divide_by_normals
  
  % N: 2 in ploidy column of sample info file (these samples may be used for normalization)
  % CTRL: control in type are normal samples used as internal control
  
  % split samples to ones used for normalization and the ones that should be normalized --> REF, CL
  ref=find(ML3.supdat(find_supid(ML3,'N'),:));
  notref=setdiff(1:size(ML3.dat,2),ref);
  
  if ignore_flagged_normals
    % Find flagged samples (those with hist_qc supdat > 0) and remove those
    % from the list of references
    hist_qc_supid=find_supid(ML3,'hist_qc');
    if ~isempty(hist_qc_supid)
      bad_ref = find(ML3.supdat(find_supid(ML3,'hist_qc'),ref)>0);
      bad_ref = ref(bad_ref);
      good_ref = setdiff(ref,bad_ref);
    else
      disp('No hist_qc label... using all reference samples');
      good_ref = ref;
      bad_ref=[];
    end
  else
    good_ref = ref;
    bad_ref=[];
  end
  
  REF=reorder_D_cols(ML3,good_ref);

  cl_idx=[]; 
  if normalize_normals
    cl_idx=[cl_idx good_ref bad_ref];
  end
  if normalize_tumors
    cl_idx=[cl_idx notref];
  end
  % CL has all samples to be normalized
  CL=reorder_D_cols(ML3,cl_idx);
  
  clear ML3
  
  % find unique CORE values
  uc=unique(CL.supdat(find_supid(CL,'CORE'),:));
  if use_all_cores_for_ref
    uc=1;
  end

  % all_center used if no normals in core
  switch norm_collapse_method
    case 'mean'
     all_center=mean(REF.dat,2);
   case 'median'
     all_center=median(REF.dat,2);
   otherwise
     all_center=mean(REF.dat,2);
  end
  
  CL.ref=zeros(size(CL.dat));
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
      % normalize using all_center
      CL.dat(:,cl_incore)=CL.dat(:,cl_incore)-repmat(all_center,1,length(cl_incore));
      CL.used_normals(cl_incore)=cellstr(repmat('center of all normals',1,length(cl_incore)));
      disp('Could not find reference for samples:');
      disp(catn(CL.sdesc(cl_incore)));
    else
      switch norm_select_method
       case 'closest_n'
        disp(['Using n_closest_normals (n=' num2str(n_closest_n) ')']);
        d=dist(REF.dat(CL.chrn~=23,ref_incore),CL.dat(CL.chrn~=23,cl_incore),'euclidean_after_med_subtract',1,1);
        [tmp,same_idx1,same_idx2]=match_string_sets_hash(REF.sdesc(ref_incore),CL.sdesc(cl_incore));
        d(sub2ind(size(d),same_idx1,same_idx2))=Inf;
        [sd,si]=sort(d,1);
        si(isinf(sd))=NaN;
        for j=1:length(cl_incore)
          if use_paired 
            % find matched normal
            use_normals=strmatch(CL.sis(cl_incore(j)).paired,REF.sdesc,'exact'); 
            if length(use_normals)>1 
              error(['more than 1 matched normal to samples: ' CL.sdesc(cl_incore(j))]);
            end
          end
          if ~use_paired || isempty(use_normals)
            use_normals=ref_incore(si(1:min(size(si,1),n_closest_n),j));
            use_normals=use_normals(~isnan(use_normals));
          end
          
          switch norm_collapse_method
           case 'mean'
            disp(['Normalizing ' CL.sdesc{cl_incore(j)} ' against mean of ' num2str(n_closest_n) ' closest normals']);
            CL.dat(:,cl_incore(j))=CL.dat(:,cl_incore(j))-mean(REF.dat(:,use_normals),2);
           case 'median'
            disp(['Normalizing ' CL.sdesc{cl_incore(j)} ' against median of ' num2str(n_closest_n) ...
                  ' closest normals']);
            
%            if strcmp(CL.sdesc{cl_incore(j)},'CHAMS_p_Sty31_(CO-124089)_Mapping250K_Sty_H02_112346')
%  who the hell put this in and why it was not taken out. (GG 2007/11/28)
%            if cl_incore(j)==3000
%              keyboard
%            end
            tmp=CL.dat(:,cl_incore(j));
            mu=median(REF.dat(:,use_normals),2);
            
            CL.dat(:,cl_incore(j))=tmp-mu; % CL.dat(:,cl_incore(j))-median(REF.dat(:,use_normals),2);
            CL.ref(:,cl_incore(j))=mu;
%            CL.dat(:,cl_incore(j))=CL.dat(:,cl_incore(j))-median(REF.dat(:,use_normals),2);
           case 'tangent'
            disp(['Using tangent normalization with ' num2str(n_closest_n) ' closest normals']);
            if length(use_normals)==1
              warning('Using tangent with 1 sample');
            end
            xx=tangent_normalization(CL.dat(:,cl_incore(j)),REF.dat(:,use_normals));
            CL.dat(:,cl_incore(j))==xx;
           otherwise 
            error('no such method');
          end
          
          CL.used_normals{j}=REF.sdesc(use_normals);
          % FIXME: DO WE NEED THIS? No!
          if (0)
            tmp=si(1:min(size(si,1),Inf),j);
            tmp(isnan(tmp))=[];
            CL.temp{j,1}=REF.sdesc(ref_incore(tmp));
            CL.temp{j,2}=sd(:,j);
          end
        end
        
       case 'all'
        tangent_ref_set=[];
        for j=1:length(cl_incore)
          if use_paired 
            % find matched normal
            if ~strmatch(lower(CL.sis(cl_incore(j)).paired),'yes')
              use_normals=strmatch(CL.sis(cl_incore(j)).paired,REF.sdesc,'exact'); 
              if length(use_normals)>1 
                error(['more than 1 matched normal to samples: ' CL.sdesc(cl_incore(j))]);
              end
            else
              use_normals=[];
            end
          end
          if ~use_paired || isempty(use_normals)
            use_normals=ref_incore;
          end
          
          use_normals=setdiff(use_normals,cl_incore(j));
          if isempty(use_normals)
            CL.dat(:,cl_incore(j))=CL.dat(:,cl_incore(j))-repmat(all_center,1,length(cl_incore));
            CL.used_normals{j}='center of all normals';
          else
            switch norm_collapse_method
             case 'mean'
              disp([ 'Core ' num2str(uc(i)) ') Dividing ' num2str(length(cl_incore)) ...
                     ' by ' num2str(length(ref_incore)) ' reference samples (mean)']);
              CL.dat(:,cl_incore(j))=CL.dat(:,cl_incore(j))-mean(REF.dat(:,use_normals),2);
             case 'median'
              disp([ 'Core ' num2str(uc(i)) ') Dividing ' num2str(length(cl_incore)) ...
                     ' by ' num2str(length(ref_incore)) ' reference samples (median)']);
              tmp=CL.dat(:,cl_incore(j));
              mu=median(REF.dat(:,use_normals),2);
              CL.dat(:,cl_incore(j))=tmp-mu; % CL.dat(:,cl_incore(j))-median(REF.dat(:,use_normals),2);
              CL.ref(:,cl_incore(j))=mu;
             case 'tangent'
              disp(['Using tangent normalization with all normals (in core)']);
              if length(use_normals)==1
                warning('Using tangent with 1 sample');
                xx=tangent_normalization(CL.dat(:,cl_incore(j)),REF.dat(:,use_normals));
                CL.dat(:,cl_incore(j))==xx;
              else
                if ~isempty(setxor(tangent_ref_set,use_normals)) % did we calc q for the same normals
                  [xx,tangent_q]=tangent_normalization(CL.dat(:,cl_incore(j)),REF.dat(:,use_normals));
                  tangent_ref_set=use_normals;
                else
                  xx=tangent_normalization(CL.dat(:,cl_incore(j)),REF.dat(:,use_normals),tangent_q);
                end
                CL.dat(:,cl_incore(j))==xx;
              end
             otherwise 
              error('no such method');
            end
            CL.used_normals{j}=REF.sdesc(use_normals);
          end
          
        end
      end
    end
  end
  if median_scale_after_norm
    disp('subtracting median of values (after dividing by normals) excluding the X');
    CL.dat=CL.dat-repmat(median(CL.dat(CL.chrn~=23,:),1),size(CL.dat,1),1);
  end
else
  CL=ML3;
  disp('Not dividing by normals');
  all_center=[];
end

CL.dat=2.^(CL.dat+1); % back to copy number



