function m=mode_freq(x,freq,dim)

if ~exist('dim','var')
  dim=1;
end

[m,f]=mode(x,dim);
f=f/size(x,dim);
m(f<freq)=NaN;


