function filter_features(in_fname,out_fname,thresh)

if ischar(thresh)
  thresh=str2num(thresh);
end

D=read_mit_gct_file(in_fname);
D1=filter_D_rows(D,struct('method','maxval','thresh',thresh));
write_mit_gct_file(out_fname,D1);

