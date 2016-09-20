function affy_info=read_affy_info(fname)

f=read_dlm_file(fname,',');
for i=2:length(f)
  j=i-1;
  affy_info(j).probeset=f{i}{1};
  affy_info(j).identifier=f{i}{2};
  affy_info(j).locus_id=f{i}{3};
  affy_info(j).name=f{i}{4};
  affy_info(j).go=f{i}{5};
  affy_info(j).prot_dom=f{i}{6};
  affy_info(j).pathway=f{i}{7};
  affy_info(j).chr=f{i}{8};
  affy_info(j).desc=f{i}{9};
end

if (0)
  u95ps=strvcat(affy_info(:).probeset);
  [us,usi]=sortrows(u95ps);
  [gs,gsi]=sortrows(gacc);
  
  max(max(abs(gs-us))) %0 -> identical
  [tmp,gsi_inv]=sort(gsi);
  u95av2_nov=u95av2(usi(gsi_inv(1:10)));
  
  save u95av2_nov.mat u95av2_nov
end
