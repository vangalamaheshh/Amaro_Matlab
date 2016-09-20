function C2=recover_peaks(C,is_log,peak_sz,min_abs_delta)

if ~exist('peak_sz','var')
  peak_sz=5;
end

if ~exist('min_abs_delta','var')
  min_abs_delta=2;
end

if mod(peak_sz,2)~=1
  error('peak_sz must be odd');
end

if (0)
  delta=C2.dat-C2.cbs;
  y=prctile(delta(:),[10 90]);
  central_dat=delta(:);
  central_dat=central_dat(find(central_dat>=y(1) & central_dat<=y(2)));
  m=mean(central_dat);
  s=std(central_dat);
  % m+10*s = 3.3189
  % m-10*s = -3.3843
  
  % we'll take +/- 3.3
  
  C2.cbs_fixed=C2.cbs;
  C2.cbs_fixed(abs(delta)>=3.3)=C2.dat(abs(delta)>=3.3);
end

C2=C;
clear C;
C2.cbs_fixed=C2.cbs;
% topn=10;

if ~is_log
  a1=C2.dat; a1(a1<0.1)=0.1; a1=log2(a1);
  a2=C2.cbs; a2(a2<0.1)=0.1; a2=log2(a2);
  delta_nonan=a1-a2;
else
  delta_nonan=C2.dat-C2.cbs;
end

%y=prctile(delta_nonan(:),[10 90]);
%central_dat=delta_nonan(:);
%central_dat=central_dat(find(central_dat>=y(1) & central_dat<=y(2)));
%m=mean(central_dat);
%s=std(central_dat);

delta_nonan(isnan(delta_nonan))=0;
X=C2;
X.dat=abs(delta_nonan);

Y=C2;
Y.dat=delta_nonan;

X1=smooth_copy_number(X,peak_sz,0,'min');% minimum of abs delta
X2=smooth_copy_number(Y,peak_sz,0,'prodrange'); % prod of min and max

for i=1:size(C2.dat,2)
  pos=find(X1.smooth(:,i)>min_abs_delta & X2.smooth(:,i)>0);
%  [ss,si]=sort(abs(delta_nonan(:,i)));
%  si=si((end-topn+1):end);
  disp([ 'recovering ' num2str([ i length(pos) ])]);
  for j=1:length(pos)
    rng=(pos(j)-(peak_sz-1)/2):(pos(j)+(peak_sz-1)/2);
    C2.cbs_fixed(rng,i)=C2.dat(rng,i);
  end
end
