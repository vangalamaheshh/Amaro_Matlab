function c=crc(c,crc_table,str)
% c=crc(c,crc_table,str)
%    calculates the CRC of the string or strings in str
%    Initialize crc_table with crc_init function.
%    c is of type uint32 and can start with 0.
%   
% Gaddy Getz
% Cancer Genomics
% The Broad Institute
% gadgetz@broad.mit.edu
%


% c=crc(c,crc_table,str)
% c is uint32. start with 0


%WORD crc_calc(WORD crc, char *buf, unsigned nbytes)
%{
%      unsigned char *p, *lim;
%      p = (unsigned char *)buf;
%      lim = p + nbytes;
%      while (p < lim)
%      {
%            crc = (crc >> 8 ) ^ crc_table[(crc & 0xFF) ^ *p++];
%      }
%      return crc;
%}

nstr=size(str,1);

if nstr>1
  if length(c)==1
    c=repmat(c,nstr,1);
  end
  
  p=uint8(str);
  ln=size(p,2);

  FF=uint32(repmat(255,nstr,1));
  for i=1:ln
    %
    %    c=bitxor(bitshift(c,-8),crc_table(bitxor(bitand(c,uint32(255)),p(i))+1));
    c1=bitshift(c,-8);
    c2=crc_table( double(bitxor( bitand(c,FF), uint32(p(:,i)) ))+1);
    c=bitxor(c1,c2);
  end  
else  
  p=uint8(str);
  ln=size(p,2);
  
  for i=1:ln
    %
    %    c=bitxor(bitshift(c,-8),crc_table(bitxor(bitand(c,uint32(255)),p(i))+1));
    c1=bitshift(c,-8);
    c2=crc_table( double(bitxor( bitand(c,uint32(255)), uint32(p(i)) ))+1);
    c=bitxor(c1,c2);
  end  
end
