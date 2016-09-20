% ---
% $Id$
% $Date$
% $LastChangedBy$
% $Rev$

function [q,p,d,ads]=snp_score_permutations(C,score_type,nperm)

if ischar(score_type)
  score_type.method=score_type;
end

switch score_type.method
 case 'nxa'
  
  if isfield(score_type,'amp_thresh')
    amp_thresh=score_type.amp_thresh;
  else  
    error('No amp_thresh');
  end
  
  if isfield(score_type,'res_per_sample')
    score_type.res=score_type.res_per_sample/size(C.dat,2);
  end
  
  if isfield(score_type,'del_thresh')
    del_thresh=score_type.del_thresh;
  else  
    error('No del_thresh');
  end
  verbose('Ignoring nperm: performing exact...',10);
 
  n=size(C.dat,2);
  s=size(C.dat,1);

  %%% Amplifications
  x=C.dat;
  x(x<amp_thresh)=0;
  if isfield(score_type,'max_segment_size')
    x=remove_large_segments(C,x,score_type.max_segment_size);
  end
  
  if nnz(isnan(x))>0
    tmp=make_D(x);
    tmp.chrn=C.chrn;
    tmp.pos=C.pos;
    tmp=impute_missing_values(tmp,'max');
    if isfield(score_type,'smooth_sz')
      tmp.cbs=tmp.dat;
      tmp=smooth_cbs(tmp,score_type.smooth_sz);
      tmp.dat=tmp.cbs;
    end
    x=tmp.dat;
  end
  
%  rla=runlength(x);
  ha=cell(1,size(x,2));
%  adsm{1}=zeros(size(x));
  for i=1:size(x,2)
%    tmp1=floor(rla{i}(:,3)/n/score_type.res)+1;
%    sz=rla{i}(:,2)-rla{i}(:,1)+1;
%    tmp2=sparse(1:length(sz),tmp1,sz);
%    ha{i}=full(sum(tmp2,1))'/s;
    ha{i}=histc(x(:,i)/n,0:score_type.res:(max(x(:,i))/n+score_type.res))/s;
    if mod(i,100)==0
      verbose(num2str(i),10);
    end
  end
  if (0)
    for i=1:size(x,2)
      ha{i}=histc(x(:,i)/n,0:score_type.res:(max(x(:,i))/n+score_type.res))/s;
      %    adsm{1}(:,i)=floor(x(:,i)/n/score_type.res);
    end
  end
  
  
  if isfield(score_type,'lsfdir')
    if isfield(score_type,'nparts')
      nparts=score_type.nparts;
    else
      nparts=10;
    end
    prts=get_parts(1:length(ha),nparts);
    hn=zeros(nparts,1);
    l=lsf(score_type.lsfdir);
    for prti=1:nparts
      [l,hn(prti)]=bsub(l,{'damp_p'},'conv_many',{ha(prts{prti})});
    end
    [l,res]=wait(l); % wait for all
    
    damp=res{hn(1)}.damp_p;
    for i=2:nparts
      damp=conv(damp,res{hn(i)}.damp_p);
    end
  else
    verbose(['Amp:.' num2str(i)],10);
    damp=ha{1};
    for i=2:length(ha);
      damp=conv(damp,ha{i});
    end
  end
  if (0)
    %%% fft
    X=zeros(max(cellfun('length',ha)),length(ha));
    for i=1:length(ha)
      X(1:length(ha{i}),i)=ha{i};
    end
    damp_fft=conv_many_fft(X);
    qq=damp-damp_fft(1:length(damp));
    close all
    plot(qq./(damp+eps))
    disp('hit a key');
    pause
  end
  
  sd=sum(damp);
  if abs(sd-1)>0.1
    error('not a distribution');
  else
    damp=damp./sum(damp);
  end
  
  % maybe replace with this
  if (0)
    m=2*round(1/score_type.res);
    F=zeros(m,n);
    for i=1:n
      F(:,i)=fft(ha{i},m);
    end
    damp2=ifft(prod(F,2));
    damp2(damp2<0)=0;
    % still not exactly the same at low values...
  end
  
  clear x;

  %%%% Deletions
  
  x=-C.dat;
  x(x<-del_thresh)=0;
  if isfield(score_type,'max_segment_size')
    x=remove_large_segments(C,x,score_type.max_segment_size);
  end  
  if nnz(isnan(x))>0
    tmp=make_D(x);
    tmp.chrn=C.chrn;
    tmp.pos=C.pos;
    tmp=impute_missing_values(tmp,'max');
    if isfield(score_type,'smooth_sz')
      tmp.cbs=tmp.dat;
      tmp=smooth_cbs(tmp,score_type.smooth_sz);
      tmp.dat=tmp.cbs;
    end
    x=tmp.dat;
  end
  
  hd=cell(1,size(x,2));
%  adsm{2}=zeros(size(x)); 
  for i=1:size(x,2)
    hd{i}=histc(x(:,i)/n,0:score_type.res:(max(x(:,i))/n+score_type.res))/s;
%    adsm{2}(:,i)=floor(x(:,i)/n/score_type.res);
    if mod(i,100)==0
      verbose(num2str(i),10);
    end
  end
  clear x;
  
  ddel=hd{1};
  verbose(['Del:.' num2str(i)],10);
  for i=2:length(ha);
    ddel=conv(ddel,hd{i});
  
  end
  
  
  if (0)
    %%% fft
    X=zeros(max(cellfun('length',hd)),length(hd));
    for i=1:length(ha)
      X(1:length(hd{i}),i)=hd{i};
    end
    ddel_fft=conv_many_fft(X);
    qq=ddel-ddel_fft(1:length(ddel));
    plot(qq./(ddel+eps))
    disp('hit a key');
    pause
  end
  
 otherwise
  error('No such score method');
end
d{1}=damp;
d{2}=ddel;
for k=1:2
  t{k}=flipud(cumsum(flipud(d{k})));
%  t{k}=1-cumsum([0; d{k}]);  % not accurate
end
[ads{1},ads{2}]=snp_score(C,score_type);
% ads2{1}=sum(adsm{1},2)+1;
% ads2{2}=sum(adsm{2},2)+1;

%for k=1:2
%  p{k}=t{k}(min(ads2{k},length(t{k})));
%  q{k}=calc_fdr_value(p{k});
%end

%if (1)
  for k=1:2
    p{k}=t{k}(min(1+floor(ads{k}/score_type.res),length(t{k})));
    q{k}=calc_fdr_value(p{k});
  end
end


