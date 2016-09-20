function D=gsea_match_sets(s,D,use_symb)

if ~exist('use_symb','var')
  use_symb=0;
end

gsupdat_T=sparse(size(D.dat,1),length(s));
for i=1:length(s)
  if mod(i,10)==0
    disp([num2str(i) '/' num2str(length(s))]);
  end
  if (i==1)
    if use_symb
      [dum,m1,m2,h,us1j]=match_string_sets_hash(D.gsymb,s(i).genes);
    else
      [dum,m1,m2,h,us1j]=match_string_sets_hash(D.gacc,s(i).genes);
    end
  else
    if use_symb
      [dum,m1,m2]=match_string_sets_hash([],s(i).genes,h,us1j);
    else
      [dum,m1,m2]=match_string_sets_hash([],s(i).genes,h,us1j);
    end
  end    
  gsupacc{i}=s(i).name;
  gsupdesc{i}=[s(i).name ' - ' s(i).type];
  if ~isempty(m1)
    [um1,um1i]=unique_keepord(m1);
    gsupdat_T(um1,i)=m2(um1i);
  else
    disp(['could not match any gene in ' s(i).name ' (#' num2str(i) ')']);
  end
end

if isfield(D,'gsupdat')
  D.gsupdat=[D.gsupdat; gsupdat_T'];
else
  D.gsupdat=gsupdat_T';
end

if isfield(D,'gsupacc')
  D.gsupacc=strvcat(D.gsupacc,strvcat(gsupacc));
else
  D.gsupacc=strvcat(gsupacc);
end

if isfield(D,'gsupdesc')
  D.gsupdesc=strvcat(D.gsupdesc,strvcat(gsupdesc));
else
  D.gsupdesc=strvcat(gsupdesc);
end
  
