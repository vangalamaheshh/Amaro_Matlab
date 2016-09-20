function D2=filter_D_rows(D,filter)

D=add_history(D,mfilename,filter);

if iscell(filter)
  D2=D;
  for i=1:length(filter)
    D2=filter_D_rows(D2,filter{i});
  end
  return
end

if ischar(filter)
  tmp.method=filter;
  filter=tmp;
end

switch filter.method
 case 'mad'
  s=mad(D.dat,1,2); % use median(abs(x-median(x)))
  filt_idx=find(s>=filter.thresh);
 case 'topmax'
  m=max(D.dat,[],2); 
  [sm,si]=sort(m,1,'descend');
  filt_idx=si(1:filter.n);
 case 'topmad'
  s=mad(D.dat,1,2); % use median(abs(x-median(x)))
  [ss,si]=sort(-s);
  filt_idx=si(1:filter.n);
 case 'std'
  s=nanstd(D.dat,0,2);
  filt_idx=find(s>=filter.thresh);
 case 'topstd'
  s=nanstd(D.dat,0,2);
  [ss,si]=sort(-s);
  filt_idx=si(1:min(filter.n,length(si)));  
 case 'toprange'
  r=range(D.dat,2);
  [rs,ri]=sort(-r);
  filt_idx=ri(1:filter.n);
 case 'top_prctile_range'
  y=prctile(D.dat',filter.prctile);
  r=diff(y);
  [rs,ri]=sort(-r);
  filt_idx=ri(1:filter.n);
 case 'prctile_range'
  y=prctile(D.dat',filter.prctile);
  r=diff(y);
  filt_idx=find(r>=filter.thresh);
 case 'minmaxval'
  m=nanmax(D.dat,[],2);
  [ss,si]=sort(m);
  filt_idx=si(1:filter.n);
 case 'minval'
  m=nanmin(D.dat,[],2);
  filt_idx=find(m>=filter.thresh);
 case 'maxval'
  m=nanmax(D.dat,[],2);
  filt_idx=find(m>=filter.thresh);
 case 'n_at_least_above'
  filt_idx=find(sum(D.dat>=filter.thresh,2)>=filter.at_least);
 case {'range','max-min'}
  filt_idx=find(max(D.dat,[],2)-min(D.dat,[],2)>=filter.thresh);
 case '2^max-2^min'
  filt_idx=find(2.^max(D.dat,[],2)-2.^min(D.dat,[],2)>=filter.thresh);
 case 'max/min'
  filt_idx=find(max(D.dat,[],2)./min(D.dat,[],2)>=filter.thresh);
 case 'strpattern'
  filt_idx=find(~cellfun('isempty',regexp(getfield(D,filter.field),filter.regexp)));
 case 'non_nan'
  filt_idx=find(~any(isnan(D.dat),2));
 otherwise
  error('no such filter');
end

if isfield(filter,'inverse') && filter.inverse==1
  filt_idx=setdiff(1:size(D.dat,1),filt_idx);
end
D2=reorder_D_rows(D,filt_idx);
D2.filt_idx=filt_idx;

