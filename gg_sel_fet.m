function [fetid,res]=gg_sel_fet(D,supid,nfet,selection_params)

res=[];
if ischar(selection_params)
  tmp.method='gp';
  tmp.score=selection_params;
  selection_params=tmp;
end

switch selection_params.method
 case 'gp'
  cls=unique_keepord(D.supdat(supid,:));
  ncls=length(cls);
  
  fetpercls=get_parts(1:nfet,ncls);
  fetid=zeros(nfet,1);
  for i=1:ncls
    global FET_SEL_SI
    %     FET_SEL_SI=[]; % FOR DEBUG: MAKE SURE FEATURES ARE CALC. EACH TIME
    if isempty(FET_SEL_SI) || ~isfield(D,'sel_fet_id') || ...
          length(FET_SEL_SI)<D.sel_fet_id || isempty(FET_SEL_SI{D.sel_fet_id}) || ...
          length(FET_SEL_SI{D.sel_fet_id})<i || isempty(FET_SEL_SI{D.sel_fet_id}{i})
      [p,s]=differential_analysis(D,find(D.supdat(supid,:)~=cls(i)), ...
                                  find(D.supdat(supid,:)==cls(i)), ...
                                  selection_params.score,1);
      [ss,si]=sort(s);
      si=flipud(si);
      fprintf(1,'x');
      if isfield(D,'sel_fet_id')
        FET_SEL_SI{D.sel_fet_id}{i}=si;
%        fprintf(1,'x');
      end
    else
      si=FET_SEL_SI{D.sel_fet_id}{i};
%      fprintf(1,'.');
    end
    fetid(fetpercls{i})=si(1:length(fetpercls{i}));
  end
 case 'gp_mult'
  cls=unique_keepord(D.supdat(supid,:));
  ncls=length(cls);
  
  fetpercls=get_parts(1:nfet,ncls);
  fetid=zeros(nfet,1);
  for i=1:ncls
    [p,s]=differential_analysis(D,find(D.supdat(supid,:)~=cls(i)),find(D.supdat(supid,:)==cls(i)),selection_params.score,1);
    [ss,si]=sort(abs(s));
    si=flipud(si);
    fetid(fetpercls{i})=si(1:length(fetpercls{i}));
  end
  
 case '2class_score'
  % assume only 2 classes
  [p,s]=differential_analysis(D,find(D.supdat(supid,:)==0),find(D.supdat(supid,:)==1),selection_params.score,1);
  [ss,si]=sort(abs(s));
  si=flipud(si);
  fetid=si(1:nfet);
  
 case '2class_pv'
  % assume only 2 classes
  global FET_SEL_SI
  %     FET_SEL_SI=[]; % FOR DEBUG: MAKE SURE FEATURES ARE CALC. EACH TIME
  if isempty(FET_SEL_SI) || ~isfield(D,'sel_fet_id') || ...
        length(FET_SEL_SI)<D.sel_fet_id || isempty(FET_SEL_SI{D.sel_fet_id}) 
    [p,s]=differential_analysis(D,find(D.supdat(supid,:)==0),find(D.supdat(supid,:)==1),selection_params.score, ...
                                0);
    [ps,si]=sort(p);
    fprintf(1,'x');
    if isfield(D,'sel_fet_id')
      FET_SEL_SI{D.sel_fet_id}=si;
    end
  else
    si=FET_SEL_SI{D.sel_fet_id};
    if isempty(si)
      keyboard
    end
    %      fprintf(1,'.');
  end
  fetid=si(1:nfet);  
  
 case '2class_pv_fdr'
  % assume only 2 classes
  global FET_SEL_SI
  %     FET_SEL_SI=[]; % FOR DEBUG: MAKE SURE FEATURES ARE CALC. EACH TIME
  if isempty(FET_SEL_SI) || ~isfield(D,'sel_fet_id') || ...
        length(FET_SEL_SI)<D.sel_fet_id || isempty(FET_SEL_SI{D.sel_fet_id}) 
    [p,s]=differential_analysis(D,find(D.supdat(supid,:)==0),find(D.supdat(supid,:)==1),selection_params.score, ...
                                0);
    fdr=calc_fdr_value(p);
    fprintf(1,'x');
    if isfield(D,'sel_fet_id')
      FET_SEL_SI{D.sel_fet_id}=fdr;
    end
  else
    fdr=FET_SEL_SI{D.sel_fet_id};
    if isempty(fdr)
      keyboard
    end
    %      fprintf(1,'.');
  end
  fetid=find(fdr<=selection_params.fdr_thresh);
 
 case 'code_entropy'

  if (0)
    X=D;  
    if isfield(selection_params,'ranges')
      X=discretize_D(X,struct('method','val','ranges',selection_params.ranges));
    end
    if isfield(selection_params,'freq')
    X=aggregate_D(X,'cols',supid,struct('method','mode_freq','freq',selection_params.freq));
    end
    X=filter_D_rows(X,'non_nan');
    filt_idx=X.filt_idx;
  else
    X=D;
    filt_idx=1:size(X.dat,1);
  end
  
  X.gorigidx=filt_idx;
  C={};
  si=[];
  for i=1:selection_params.n_codes
    disp(['CODE #' num2str(i)]);
    C{i}=make_cipher(X,selection_params.use_base_entropy);
    C{i}=add_code_dend(C{i});
    si=[si; as_column(C{i}.D.gorigidx)];
    keepi=setdiff(1:size(X.dat,1),C{i}.idx);
    X=reorder_D_rows(X,keepi);
  end
  fetid=si;
  res=C;

 otherwise
  error('No such method');
end

% disp(['selected ' num2str(fetid') ' from  1:' num2str(length(ss)) ]);


