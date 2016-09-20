function D=collapse_supdat(D,acc,desc,range,remove_flag)

 s=sum(D.supdat(range,:),1);
 if any(s>1)
   disp('Not unique');
 else
   if ~isempty(find(s==0))
     [tmp,ind]=max([ 1-sum(D.supdat(range,:)); D.supdat(range,:)]);
   else
    [tmp,ind]=max([D.supdat(range,:)]);
  end
  D.supdat(end+1,:)=ind;
  if ~isempty(find(s==0))  
    longname=[ desc ': 1-Other'];
    longacc=[ acc ': 1-OTHER'];
    for i=1:length(range)
      longacc=[longacc '/' num2str(i+1) '-' deblank(D.supacc(range(i),:))];
      longname=[longname '/' num2str(i+1) '-' deblank(D.supdesc(range(i),:))];
    end
  else
    longname=[ desc ': ' '1-' deblank(D.supdesc(range(1),:))];
    longacc=[ acc ': ' '1-'  deblank(D.supacc(range(1),:))];
    for i=2:length(range)
      longacc=[longacc '/' num2str(i) '-' deblank(D.supacc(range(i),:))];
      longname=[longname '/' num2str(i) '-' deblank(D.supdesc(range(i),:))];
    end
  end
  D.supdesc=strvcat(D.supdesc,longname);
  D.supacc=strvcat(D.supacc,longacc);
  if (remove_flag)
    D.supdat(range,:)=[];
    D.supdesc(range,:)=[];
    D.supacc(range,:)=[];
  end
end

sdm=mat2cell(D.supdat,ones(size(D.supdat,1),1),ones(size(D.supdat, ...
                                                  2),1));
for i=1:size(D.supacc)
  pos=find(D.supacc(i,:)==':');
  if ~isempty(pos)
    sacc{i}=deblank(D.supacc(i,1:(pos(1)-1)));
  else
    sacc{i}=deblank(D.supacc(i,:));    
  end
end

D.sup=cell2struct(sdm,sacc,1);
