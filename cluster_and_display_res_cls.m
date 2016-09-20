function cluster_and_display_res_cls(resfile,clsfile,sparams, ...
                                     gparams,disptype,colortab,outfname,printparams)

if strcmp(file_ext(resfile),'gct')
  D=read_mit_gct_file(resfile);
else
  D=read_mit_res_file(resfile);
end

D=read_mit_new_cls_file(D,clsfile);
D=add_supmark(D,colortab);
[Dord,gdend,sdend]=two_way_clustering(D,gparams,sparams);

display_D(Dord,gdend,sdend,disptype);
print_D(outfname,printparams,1,[]);
