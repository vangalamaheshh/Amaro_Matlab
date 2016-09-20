function idx=findstrings_list(strs,strlist)
% find indices of strlist in strs
% returns empty cell if not found
% assumes strlist and strs are unique

if(0)
  nstrs=size(strs,1);
  [sorted_strs,si]=sortrows(strs);
  [tmp,rev_si]=sort(si);
  
  [ustrs,ui,uj]=unique(strvcat(sorted_strs,strlist),'rows');
  
  if size(ustrs,1)~=nstrs
    disp([ 'strlist contains stings that dont match any in strs or '...
           'strs is not unique']);
  end

  idx=si(uj((nstrs+1):end));
else
  ns1=size(strs,1);
  [us1,us1i,us1j]=unique(strs,'rows');
  nus1=size(us1,1);
  
  ns2=size(strlist,1);
  [us2,us2i,us2j]=unique(strlist,'rows');
  nus2=size(us2,1);
  
  [ustrs,ui,uj]=unique(strvcat(us1,us2),'rows');
  uj1=uj(1:nus1);
  uj2=uj(nus1+(1:nus2));

  U1=cell(nus1,1);
  U2=cell(nus2,1);
  if nus1<nus2
    for i=1:nus1
      tmp=find(uj2==uj1(i));
      U1{i}=tmp;
      if ~isempty(tmp);
        U2{tmp}=i;
      end
    end
  else
    for i=1:nus2
      tmp=find(uj1==uj2(i));
      U2{i}=tmp;
      if ~isempty(tmp);
        U1{tmp}=i;
      end
    end
  end
  
  if (ns1 ~= nus1) | (ns2 ~= nus2)
    error(['either strs or strlist are not unique']);
  else
    idx=U2(us2j);
    nonempty=find(~cellfun('isempty',idx));
    tmp=us1i(cell2mat(idx(nonempty)));
    for i=1:length(nonempty)
      idx{nonempty(i)}=tmp(i);
    end
  end 
end




