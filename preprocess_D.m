function D=preprocess_D(D,params)

D=add_history(D,mfilename,params);

if ~exist('params','var')
  params='';
end

if ischar(params)
  params1.method=params;
  params=params1;
end

if ~isfield(params,'dont_save_before_preproc') || ~params.dont_save_before_preproc
  D.before_preproc_dat=D.dat;
end

if iscell(params)
  for i=1:length(params)
    D=preprocess_D(D,params{i});
  end
else
  switch params.method
   case 'confounded_preproc'
    u=unique_keepord(D.supdat(params.supid,:));
    origidx=[];
    for i=1:length(u)
      ci=find(D.supdat(params.supid,:)==u(i));
      origidx=[origidx ci];
      Ds{i}=preprocess_D(reorder_D_cols(D,ci),params.preproc_params);
    end
    D2=unite_Ds(Ds,'cols');
    %  D2=rmfield(D2,'origidx');
    %  D2=rmfield(D2,'residx');
    [dum,revord]=sort(origidx);
    D=reorder_D_cols(D2,revord);
   case 'row_center_and_normalize'
    if isfield(params,'center') & isfield(params,'scale')
      D.dat=D.dat-repmat(params.center,1,size(D.dat,2));
      D.dat=D.dat./repmat(params.scale,1,size(D.dat,2));
    else
      [D.dat,D.preproc_center,D.preproc_scale]= ...
          dna_norm(D.before_preproc_dat);
    end
   case 'row_median_and_mad'
    if isfield(params,'median') & isfield(params,'mad')
      D.dat=D.dat-repmat(params.median,1,size(D.dat,2));
      D.dat=D.dat./repmat(params.mad,1,size(D.dat,2));
    else
      med=median(D.before_preproc_dat,2);
      s=mad(D.before_preproc_dat,1,2);
      D.median=med;
      D.mad=s;
      D.dat=(D.dat-repmat(med,1,size(D.dat,2)))./repmat(s,1,size(D.dat,2));
    end
    
   case 'row_center'
    if isfield(params,'center')
      D.dat=D.dat-repmat(params.center,1,size(D.dat,2));
    else
      D.preproc_center=mean(D.dat,2);
      D.dat=D.dat-repmat(D.preproc_center,1,size(D.dat,2));
    end
   case 'row_center_and_normalize_unit'
    [D.dat,D.preproc_center,D.preproc_scale]=dna_norm(D.before_preproc_dat);
    sc=sqrt(size(D.before_preproc_dat,2)-1);
    D.dat=D.dat/sc;
    D.preproc_scale=D.preproc_scale*sc;
   case 'col_center_and_normalize'
    D.dat=(dna_norm(D.before_preproc_dat'))';
   case 'col_center'
    D.dat=D.before_preproc_dat-repmat(nanmean(D.before_preproc_dat,1),size(D.before_preproc_dat,1),1);
   case 'col_center_and_normalize_unit'
    D.dat=dna_norm(D.before_preproc_dat')'./sqrt(size(D.before_preproc_dat,1)-1);
   case 'rowcolumn'
    if ~ischar(params) && strcmp(params.start_with,'rows')
      D.dat=row_and_column_norm(D.before_preproc_dat); % start with rows
    else
      D.dat=row_and_column_norm(D.before_preproc_dat')'; % start with column
    end 
   case 'ceiling'
    D.dat(D.dat>=params.ceiling)=params.ceiling;
   case 'floor'
    D.dat(D.dat<=params.floor)=params.floor;
   case 'row_median_center'
    if ~isfield(params,'min_sz') || (isfield(params,'min_sz') && size(D.dat,2)>=params.min_sz)
      disp(size(D.dat,2));
      D.dat=D.dat-repmat(median(D.dat,2),1,size(D.dat,2));
    end
   case 'log_transform'
    D=log_transform(D);
   case 'rank_transform'
    [tmp,tmpi]=sort(D.dat);
    [tmp,D.dat]=sort(tmpi);
   case 'asinh'
    D.dat=asinh(D.dat/2)/log(2);
   case 'fet_scale'
    if isfield(D,'fet_scale_min')
      warning('Fetures are already scaled');
    end
    if isfield(params,'min') && isfield(params,'max')
      D.fet_scale_min=params.min;
      D.fet_scale_max=params.max;
    else
      D.fet_scale_min=min(D.dat,[],2);
      D.fet_scale_max=max(D.dat,[],2);
    end
    D.fet_scale_range=params.range;
    D.dat=(D.dat-repmat(D.fet_scale_min,1,size(D.dat,2)))./repmat(D.fet_scale_max-D.fet_scale_min,1,...
                                                      size(D.dat,2));
    D.dat=D.dat*(D.fet_scale_range(2)-D.fet_scale_range(1))+D.fet_scale_range(1);
   case 'divide_by_pivot'
    if isfield(params,'pivot')
      pivot=params.pivot;
    else
      pivot=make_pivot(D,params.pivot_params);
    end
    D.pivot=pivot;
    D.dat=D.dat-repmat(pivot.dat,1,size(D.dat,2));
   case 'divide_by_pivot_mfc'
    if isfield(params,'pivot')
      pivot=params.pivot;
    else
      [pivot,pivot_idx]=make_pivot(D,params.pivot_params);
    end
    pivot=preprocess_D(pivot,{struct('method','floor','floor',20),...
                        struct('method','ceiling','ceiling',16000)});
    if isfield(params,'remove_pivots') && params.remove_pivots
      D=reorder_D_cols(D,setdiff(1:size(D.dat,2),pivot_idx));
    end
    D=preprocess_D(D,{struct('method','floor','floor',20),...
                      struct('method','ceiling','ceiling',16000)});
    D.pivot=pivot;
    u=unique(D.supdat(params.supid,:));
    dat=repmat(D.dat,[1 1 size(pivot.dat,2)]);
    dat=dat./repmat(permute(pivot.dat,[1 3 2]),[1 size(dat,2) 1]);
    dat=log(dat);
    if (0)
      for i=1:length(u)
        uidx=find(D.supdat(params.supid,:)==u(i));
        j=find(pivot.supdat==u(i));
        dat(:,uidx,j)=NaN;
      end
    end
    dat=nanmean(dat,3);
    D.dat=dat;
   case 'threshold'
    D=threshold(D,params.val);
   case 'median_polish'
    
   case {'','none'}
    
   otherwise
    error(['no preprocessing method:' params.method]);
  end
end


