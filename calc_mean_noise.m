function [my,sy]=calc_mean_noise(mLU,sLU,hx,hdelta)

for i=1:length(hx)
  idx=find((mLU>hx(i)-hdelta) & (mLU<hx(i)+hdelta));
  my(i)=mean(sLU(idx));
  sy(i)=std(sLU(idx));
end
