function supdat2cls(D,stypes,fprefix)

for i=1:length(stypes)
  D2=D;
  si=stypes(i);
  ntps=unique(D.supdat(si,~isnan(D.supdat(si,:))));
  if length(ntps)<2  % 0 or 1
     disp(['Not enough diversity in type ' num2str(si)]);
  else
    if ~is_multiple_types(D,si) & (length(ntps)==2) 
      fname=[fprefix '.' deblank(D.supacc(si,:)) '.cls' ];
      D2.supdat=[ D2.supdat(si,:); 1-D2.supdat(si,:)];
      D2.supacc=strvcat(deblank(D2.supacc(si,:)),['not-' deblank(D2.supacc(si,:))]);
      D2.supdesc=strvcat(deblank(D2.supdesc(si,:)),['not-' deblank(D2.supdesc(si,: ...
                                                        ))]);
      if D2.supdat(2,1)==1
        D2.supdat=flipud(D2.supdat);
        D2.supdesc=flipud(D2.supdesc);
        D2.supacc=flipud(D2.supacc);
      end
%      keyboard
      verbose(['Writing to file ' fname],2);
      write_mit_cls_file(fname,D2,[1 2]);  
    else % multiple types
      [ typeacc, typedesc, D2, rg ]=decollapse_supdat(D,si);
      fname=[fprefix '.' typeacc '.cls' ];
      verbose(['Writing to file ' fname],2);
      write_mit_cls_file(fname,D2,rg);        
    end
  end
end

