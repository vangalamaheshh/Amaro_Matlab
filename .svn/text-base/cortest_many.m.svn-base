function [p,ratio] = cortest_many(dat,vec,tail,nperm)

vec=dna_norm(vec)./sqrt(size(vec,2)-1);
dat=dna_norm(dat)./sqrt(size(dat,2)-1);

C=dat*vec';
n=length(vec);
for i=1:nperm
  rv=vec(randperm(n));
  Cr(:,i)=dat*rv';
end



