function crc_table=init_crc
% c=init_crc
%    Calculates a CRC table. Used in crc function.
%   
% Gaddy Getz
% Cancer Genomics
% The Broad Institute
% gadgetz@broad.mit.edu
%

crc_table=uint32(zeros(256,1));

bufsiz=(16*1024);

for i=0:255
  k=uint32(hex2dec('C0C1'));
  j=1;
  for j1=1:8
    if bitand(i,j)
      crc_table(i+1)=bitxor(crc_table(i+1),k);
    end
    k = bitxor(bitshift(k,1),hex2dec('4003'));
    j=j*2;
  end
end

