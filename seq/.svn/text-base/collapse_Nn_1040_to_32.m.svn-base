function Nn = collapse_Nn_1040_to_32(in1,in2)
% rows:  1040 categories as in /xchip/tcga_scratch/lawrence/db/allcateg/categs.txt
% columns:  [N ->A ->C ->G ->T]
% Mike Lawrence 2009-12-11

if ~exist('in2','var')
  Nn = in1;
else
  Nn = [in1 in2];
end

if size(Nn,1)~=1040, error('input must have 10404 rows'); end
if size(Nn,2)~=5, error('input must have 5 columns (N A C G T)'); end

Nn = collapse_1040_to_128(Nn);
Nn = Nn(1:64,:)+Nn(65:128,:);
Nn = collapse_Nn_64_by_strand(Nn);

