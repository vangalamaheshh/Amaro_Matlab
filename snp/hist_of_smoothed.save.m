function [h,s,peaks,C1]=hist_of_smoothed(C,idx,rng,alpha,ignore_peak_level)

if ~exist('ignore_peak_level','var')
  ignore_peak_level=1e-6;
end

if ~exist('idx','var') || isempty(idx)
  idx=1:size(C.dat,2);
end

C=reorder_D_cols(C,idx);
if mean(C.raw(:,1))>1
  disp('taking log of C.raw');
  C.raw(C.raw<0.1)=0.1;
  C.raw=log2(C.raw)-1;
end

C=reorder_D_rows(C,C.chrn~=23);
disp('removing chrX');

if isfield(C,'medians')
  med=C.medians;
  tmp=derunlength(C.cbs_rl);
  med2=median(tmp(1:size(C.dat,1)));
  if any(med~=med2)
    error('medians do not match');
  end
else
  disp('calculating medians');
  tmp=derunlength(C.cbs_rl);
  med=median(tmp);
  C.medians=med;
end

disp('subtracting median from raw data');
r=runlength(C.dat,C.chrn);
C.raw=C.raw-repmat(C.medians,size(C.raw,1),1);
rng=[ rng Inf];


for si=1:length(idx)
  s=zeros(size(r,1),2);
  h=zeros(length(rng)-1,size(r,1));
  for i=1:size(r,1)
    s(i,1)=median(C.raw(r(i,1):r(i,2),1),1);
    s(i,2)=mad(C.raw(r(i,1):r(i,2),1),1)*mad_factor/sqrt(r(i,2)-r(i,1)+1)+alpha;
    h(:,i)=normcdf(rng(2:end),s(i,1),s(i,2))'-normcdf(rng(1:(end-1)),s(i,1),s(i,2))';
    %  h(:,i)=h(:,i)*(r(i,2)-r(i,1)+1);
  end
  hst=sum(h,2);
  hst=hst./sum(hst);
end
%clf;
%subplot(1,2,1);

spikes=histc(C.dat(:,1),rng);
plot(rng,spikes/max(spikes)*max(hst),'-k'); hold on

midrng=(rng(1:(end-1))+rng(2:end))/2;
plot(midrng,hst); hold on

%subplot(1,2,2);
%imagesc(h);

% keyboard
hst=sum(h,2);
hst=hst./sum(hst);
peaks=find_peaks(hst,ignore_peak_level);

if ~isempty(peaks)
  peaks=midrng(peaks);
  ax=axis;
  for i=1:length(peaks)
    line([peaks(i) peaks(i)],ax(3:4),'Color','r');
  end
  C1=C;
  r1=r;
  for i=1:size(r1,1)
    dst=((peaks-s(i,1))./s(i,2)).^2;
    [mdst,mdi]=min(dst);
    r1(i,3)=mdi;
  end
  tmp=derunlength(r1);
  for i=1:length(peaks)
    tmpi=find(tmp==i);
    C1.dat(tmpi)=median(C1.raw(tmpi,1));
  end
  C1.level=tmp;
end




