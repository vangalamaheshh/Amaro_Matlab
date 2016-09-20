function [D,supids]=break_marker_to_gsupdat(D,chip_type)

if ischar(chip_type)
  tmp.method=chip_type;
  chip_type=tmp;
end
supids=[];

switch chip_type.method
 case 'SNP60'
  mst=strvcat(D.marker);
  v=nan(1,size(D.dat,1));
  v(mst(:,1)=='S')=1;
  v(mst(:,1)=='C')=2;
  [D,sid]=add_D_sup(D,'MARKER: 1-CN/2-SNP','MARKER: 1-CN/2-SNP',v,'row',1);
  supids=[supids sid];
  
  n=str2int_matrix(mst);
  [D,sid]=add_D_sup(D,'MID','MID',n','row',1);
  supids=[supids sid];
 otherwise
  error('no such method');
end
