function [p,C] = corrtest_many(dat,vec,tail,nperm)

vec=dna_norm(vec)./sqrt(size(vec,2)-1);
dat=dna_norm(dat)./sqrt(size(dat,2)-1);

C=dat*vec';

if nargin>3
  n=size(vec,2);
  Cr=zeros(size(dat,1),size(vec,1));
  M=zeros(size(dat,1),size(vec,1));
  for i=1:nperm
    rv=vec(:,randperm(n));
    Cr=dat*rv';
    if mod(i,100)==0
      disp(i)
    end
    if tail==0
      M=M+(abs(Cr)>=abs(C));
    elseif tail==1
      M=M+(Cr>=C);
    end
  end
  
  p=M/nperm;
end
