
function snp_idx = snp_cnv_idx(C,cnvs_file)
  

  cnvs = read_dlm_file(cnvs_file);
  % cnvs=cat(1,cnvs{1})

  for i = 1:length(cnvs)
    for j = 1:3
      cnvs1(i,j) = str2num(char(cnvs{i}(j)));
    end
  end

  snp_idx = [];
  for i = 1:length(cnvs1)
    tmp=find(cnvs1(i,1)==C.chrn & cnvs1(i,2)<=C.pos & cnvs1(i,3)>=C.pos);
    snp_idx = cat(1,snp_idx,tmp);
  end

