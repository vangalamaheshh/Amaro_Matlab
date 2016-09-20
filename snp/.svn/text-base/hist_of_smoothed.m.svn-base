function [C,H]=hist_of_smoothed(C,rng,alpha,ignore_peak_level)

if ~exist('ignore_peak_level','var')
  ignore_peak_level=1e-6;
end

rng=[ rng rng(end)+(rng(end)-rng(end-1))]; % FIXME: remove first and last
C.hists=zeros(length(rng)-1,size(C.dat,2));
C.spikes=zeros(length(rng)-1,size(C.dat,2));
C.hist_rng=rng;
H=cell(size(C.dat,2),1);
for si=1:size(C.dat,2)
  disp(si);
  r=C.cbs_rl{si};
  s=C.rl_std{si}+alpha; % increase standard error by alpha
  h=zeros(length(rng)-1,size(r,1));
  for i=1:size(r,1)
    h(:,i)=normcdf(rng(2:end),r(i,3),s(i))'-normcdf(rng(1:(end-1)),r(i,3),s(i))';
    %  h(:,i)=h(:,i)*(r(i,2)-r(i,1)+1);
  end
  hst=sum(h,2);
  hst=hst./sum(hst);
  H{si}=h;
  C.hists(:,si)=hst;
  midrng=(rng(1:(end-1))+rng(2:end))/2;

  tmp=histc(C.dat(:,si),midrng);
  C.spikes(:,si)=tmp;

  peaks=find_peaks(hst,ignore_peak_level);
  peaks=midrng(peaks);
  C.peaks{si}=peaks;

  r1=r;
  for i=1:size(r,1)
    dst=((peaks-r(i,3))./s(i)).^2;
    [mdst,mdi]=min(dst);
    r1(i,3)=mdi;
  end
  tmp=derunlength(r1);
  for i=1:length(peaks)
    in_peak=find(tmp==i);
    C.dat(in_peak,si)=median(C.raw(in_peak,si));
  end
  C.level(:,si)=tmp;
end





