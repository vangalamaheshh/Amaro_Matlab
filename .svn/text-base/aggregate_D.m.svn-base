function D=aggregate_D(D,rc,supid,aggtype)

if ischar(aggtype)
  tmp.method=aggtype;
  aggtype=tmp;
end

if strcmp(lower(aggtype.method),'none')
  return
end

remove_supid=[];
if length(supid)>1
  v=supid;
  [D,supid]=add_D_sup(D,'TEMP','TEMP',v,rc);
  remove_supid=supid;
else
  if is_col(rc)
    v=D.supdat(supid,:);
  else
    v=D.gsupdat(supid,:); %assumes this struct of gsupdat (FIXME)
  end
end

D=reorder_D(D,rc,find(~isnan(v)));
v=v(~isnan(v));
[sv,svidx]=sort(v);
D=reorder_D(D,rc,svidx);

[dum,revidx]=sort(svidx);

[u,ui,uj]=unique(sv);

if is_col(rc)
  x=D.dat;
else
  x=D.dat';
end
assert(size(x,2)==length(sv),'sizes do not match');

y=cell(length(u),1);
for i=1:length(u)
  y{i}=x(:,find(uj==i));
end

switch aggtype.method
 case 'sum'
  yy=cellfun_any('sum(x,2)',y);
 case 'mean'
  yy=cellfun_any('mean(x,2)',y);
 case 'median'
  yy=cellfun_any('median(x,2)',y);
 case 'mode_freq'
  yy=cellfun_any(['mode_freq(x,' num2str(aggtype.freq) ',2)'],y);
 otherwise
  error('no such method');
end


if isfield(aggtype,'keepsize') && aggtype.keepsize==1
  for i=1:length(yy)
    yy{i}=repmat(yy{i},1,length(find(uj==i)));
  end
else % aggregate supdat or gsupdat
  if is_col(rc)
    if isfield(D,'supdat')
      sy=cell(length(u),1);
      for i=1:length(u)
        sy{i}=D.supdat(:,find(uj==i));
      end
      syy=cellfun_any('round(mean(x,2))',sy);      
      D.supdat=cat(2,syy{:});
    end
  else
    if isfield(D,'gsupdat')
      sy=cell(length(u),1);
      for i=1:length(u)
        sy{i}=D.gsupdat(:,find(uj==i));
      end
      syy=cellfun_any('round(mean(x,2))',sy);      
      D.gsupdat=cat(2,syy{:});
    end
  end
end
yy=cat(2,yy{:});

if is_col(rc)
  D.dat=yy;
else
  D.dat=yy';
end

if isfield(aggtype,'keepsize') && aggtype.keepsize==1
  D=reorder_D_cols(D,revidx);
elseif is_col(rc)
  [typeacc,typedesc,Dtmp,range,non_empty]=decollapse_supdat(D,supid,D.supdat(supid,:));
  D.sdesc=cellstr(Dtmp.supdesc(range,:));
else
  disp('FIX ME: add new names for rows!!!');
end

if ~isempty(remove_supid)
  if is_col(rc)
    D=reorder_D_sup(D,'col',setdiff(1:size(D.dat,2),remove_supid));
  else
    D=reorder_D_sup(D,'row',setdiff(1:size(D.dat,1),remove_supid));
  end
end

% FIX ME:
% generate new names, sup for the aggregated D
