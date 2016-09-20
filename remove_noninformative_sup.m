function D=remove_noninformative_sup(D)

v=version;
if v(5)=='0' % 7.0.0
  x=nanstd(D.supdat')';
  y=find( nansum( (abs(D.supdat-round(D.supdat))>0)' )==0);
  y=y';
else
  x=nanstd(D.supdat,[],2);
  y=find(nansum(abs(D.supdat-round(D.supdat))>0,2)==0); 
end
% remove non-integer values
D=reorder_D_sup(D,'col',intersect(find((~isnan(x) & x>0)),y));


if (0) % keep code for later use
  for i=1:size(D.supdat,1)
    notnan=find(~isnan(D.supdat(i,:)));
    u=unique(D.supdat(i,notnan));
    [st,sn]=break_sup_names(deblank(D.supacc(i,:)));
    if length(u)~=length(sn)
      ok=0;
      if ~isempty(find(u==0)) & length(u)==length(sn)+1
        ok=1;
      end
      if ~ok
        disp(['no match in ' num2str(i) ' ' st]);
      end
    end
  end
end
