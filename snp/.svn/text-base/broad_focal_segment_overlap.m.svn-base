function T = broad_focal_segment_overlap(D,Q,cyto,focal_cutoff,broad_cutoff,thresh_focal, thresh_broad,is_amp,norm_type,ref_length)

  %% Check inputs
  varlist1 = {'D','Q','cyto','focal_cutoff','broad_cutoff', ...
              'thresh_focal','thresh_broad','is_amp','norm_type','ref_length'};
  defaults = {'ERR','ERR','ERR','0.99','0.99','0.1','0.1','1','1','2'};
  required = [1,1,1,0,0,0,0,0,0,0];
  
  for idx = 1:length(varlist1)
    if required(idx) && (~exist(varlist1{idx},'var') || eval(['isempty(' varlist1{idx} ')']))
      error('Required input %s undefined.',varlist1{idx})
    elseif ~exist(varlist1{idx},'var') || eval(['isempty(' varlist1{idx} ')'])
      eval([varlist1{idx} '=' defaults{idx}]);
    end
  end
  
  disp('Initializing variables...')
  D.dat = single(D.dat*(2*is_amp-1));
  focals = zeros(size(D.dat,1),size(D.dat,2));
  focals = single(focals);
  broads = focals;
  unamp = focals;
  deleted = focals;

  disp('Normalizing by arm length...') 
  fract_chr_arm = normalize_by_arm_length(D,Q,cyto,norm_type,ref_length);
  
  focal_segs = intersect(find(fract_chr_arm < focal_cutoff),find(Q(:,4)>thresh_focal));
  broad_segs = intersect(find(fract_chr_arm > broad_cutoff),find(Q(:,4)>thresh_broad));
  
  for i=1:length(focal_segs)
    focals(Q(focal_segs(i),2):Q(focal_segs(i),3),Q(focal_segs(i),5)) = 1;
  end
  
  for i=1:length(broad_segs)
     broads(Q(broad_segs(i),2):Q(broad_segs(i),3),Q(broad_segs(i),5)) = 1;
  end
  
  figure(1)
  imagesc(focals); colorbar
  figure(2)
  imagesc(broads); colorbar
  
  no_broads = find(sum(broads,1) == 0);
  with_broads = setdiff(1:size(broads,2),no_broads);
  unamp(find(abs(D.dat)<thresh_broad))=1;
  deleted(find(D.dat<(-1*thresh_broad)))=1;
  
  good_broads = broads(:,with_broads); clear broads;
  good_focals = focals(:,with_broads); clear focals
  good_unamp=unamp(:,with_broads); clear unamp
  good_deleted=deleted(:,with_broads); clear deleted
  
  disp('Computing overlap...')
  T = zeros(2);
  
  T(1,1) = sum(sum(good_broads.*good_focals)); 
  
  T(1,2) = sum(sum(good_broads.*(ones-good_focals)));
  
  T(2,1) = sum(sum((ones-good_broads).*good_focals));
  
  T(2,2) = sum(sum(good_unamp));
  
%  % Sanity check
%  total_snps=size(good_broads,1)*size(good_broads,2);
%  
%  del_snps = sum(sum(good_deleted));
%  
%  ttotal_snps = botleft+botright+topleft+topright+del_snps;
