function D=binarize_D(D,bin_type)

if ischar(bin_type)
  bin_type.method=bin_type;
end

switch bin_type.method
 case 'threshold'
  D.dat=(D.dat>bin_type.thresh);
 otherwise
  error('no such method');
end

