function CN=read_copy_number_file(fname)

[f,fid]=read_dlm_file(fname,char(9),2);
if isempty(f{2}{end})
   tab=1;
else
   tab=0;
end 
form=[ '%s%s' repmat('%f',1,length(f{2})-2)];
CN_dat=textscan(fid,form);
CN.marker=CN_dat{1};
CN.chr=CN_dat{2};
CN.pos=CN_dat{3};
CN.cM=CN_dat{4};
CN.score=CN_dat{5};
CN.dat=cat(2,CN_dat{6:(end-tab)}); 
CN.sdesc=f{2}(6:(end-tab));
fclose(fid);
       
 
