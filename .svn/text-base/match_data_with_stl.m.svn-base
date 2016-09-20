function [D,ord,msort]=match_data_with_stl(D,stl,field)

matchvec=zeros(size(D.sdesc,1),1);
use_idx=min(findstrings({field.abr},'USE'));
if isempty(use_idx)
  error('Did not find USE column');
end
for i=1:size(stl,1)
  if stl{i,use_idx}>0 % use
    pos=strmatch(stl{i,1},D.sdesc,'exact');
    if isempty(pos)
      verbose(['*********  no match for ' stl{i,1}]);
    elseif length(pos)==1
      matchvec(pos)=i;
    else 
      verbose(['more than one match for ' stl{i,1} ' : ' ...
               num2str(pos')]);
    end
  else
    verbose(['not using ' stl{i,1}]);
  end
end


fm=find(matchvec>0);
m=matchvec(fm);
[msort,ord]=sort(m);
ord=fm(ord);
D=reorder_D_cols(D,ord);

for i=1:length(field)
  if isempty(field(i).abr) 
    break;
  end
end
i=i-1;
if i>=use_idx
  D.supacc=strvcat(field(use_idx:i).abr);
  D.supdesc=strvcat(field(use_idx:i).full);
  D.supdat=cell2mat(stl(msort,use_idx:i))';
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

D.sdesc=strvcat(stl{msort,3});

D.scans=cell(size(D.sdesc,1),1);

for j=1:2
  sjchip_field=findstrings(strvcat_empty(field(:).abr),['S' num2str(j) 'CHIP']); assert(length(sjchip_field)<=1);
  sjname_field=findstrings(strvcat_empty(field(:).abr),['S' num2str(j) 'NAME']); assert(length(sjname_field)<=1);
  for i=1:length(msort)
    if ~isempty(sjchip_field) && ~isempty(sjname_field)
      sjc=stl{msort(i),sjchip_field};
      sjn=stl{msort(i),sjname_field};
      if ~isempty(sjn) && ismember(sjn(1),'A':'Z') %scans must have first capital letter 
        D.scans{i}{j}.name=sjn;
        D.scans{i}{j}.chip=sjc;
      end
    end
  end
end


gcm_name=findstrings(strvcat_empty(field(:).abr),'GCMNAME');
if ~isempty(gcm_name)
  for i=1:length(msort)
    D.gcm_name{i}=stl{msort(i),gcm_name};
  end
end
