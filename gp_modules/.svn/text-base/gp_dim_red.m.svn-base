function gp_dim_red(in_gct,out_gct,dim_red_type,k,rc)

if ischar(dim_red_type) && ~isempty(str2num(dim_red_type))
  dim_red_type=str2num(dim_red_type);
end

if isnumeric(dim_red_type)
  params.method=enum_param(dim_red_type,...
                           {0,''; 1,'PCA';2,'SVD'; 3,'MDS'});
else
  params.method=dim_red_type;
end

if ischar(rc) && ~isempty(str2num(rc))
  rc=str2num(rc);
end

if isnumeric(rc)  
    rc=enum_param(rc,{0,'row';1,'col'});
end

if ischar(k)
  params.k=str2num(k);
end

% read gctfile
D=read_mit_gct_file(in_gct);

R=reduce_dim_D(D,rc,params.k,params);

write_mit_gct_file(out_gct,R);
