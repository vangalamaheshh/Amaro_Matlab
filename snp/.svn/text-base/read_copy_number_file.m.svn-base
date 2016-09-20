function CN=read_copy_number_file(fname,header_lines,is_loh)

if ~exist('header_lines','var')
  header_lines=1;
end
if ~exist('is_loh','var')
  is_loh=0;
end

[f,fid]=read_dlm_file(fname,char(9),1+header_lines);
if isempty(f{1+header_lines}{end})
   no_tab=1
else
   no_tab=0;
end 
if is_loh
  form=[ '%s%s' repmat('%s',1,length(f{1+header_lines})-2+no_tab)];
else
  form=[ '%s%s' repmat('%f',1,length(f{1+header_lines})-2+no_tab)];
end
CN_dat=textscan(fid,form,'bufSize',10000000,'treatAsEmpty',{'NA'});
CN.marker=CN_dat{1};
CN.chr=CN_dat{2};
CN.pos=CN_dat{3};
CN.cM=CN_dat{4};
CN.score=CN_dat{5};
keyboard
CN.dat=cat(2,CN_dat{6:(end-1)}); 
CN.sdesc=f{1+header_lines}(6:(end-1));
fclose(fid);
       
 
