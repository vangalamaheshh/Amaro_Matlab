function [s,e,sample_s,sample_e]=add_wide_peak_k(z,chr_zero,k)

rl=runlength(z);
zsc=nanmean(z,2);
[mx,mi]=max(zsc);
rg=find(zsc==mx);
st=min(rg);
en=max(rg);
samples=find(z(mi,:));

rli=find_rl(rl,mi);

if nnz(isnan(z))>0
  error('has nan');
end

zorig=z;
zi=[1; 1+find(any(diff(zorig,1,1),2))];
zj=cumsum([1; any(diff(zorig,1,1),2)]);

z=z(zi,:);

zc=z(:,setdiff(1:size(z,2),samples));
zv=z(:,samples);
zsc_c=sum(zc,2);
zsctot=zsc_c+sum(zv,2);

[mx,mi]=max(zsctot);
rg=find(zsctot==mx);
sts=min(rg);
ens=max(rg);

j=[];
if k>0
  if k==1
    [s1,e1,sample_s,sample_e]=add_wide_peak(zv,0,1);
    sts(end+1)=s1;
    ens(end+1)=e1;
  elseif k==2
    [s1,e1,smp_s1,smp_e1]=add_wide_peak(zv,0,1);
    idx_s1=setdiff(1:size(zv,2),smp_s1);
    if ~isempty(idx_s1)
      zv_s1=zv(:,idx_s1);
      [s2s,e2s,smp_s2s,smp_e2s]=add_wide_peak(zv_s1,0,1);    
      sample_s=[smp_s1 idx_s1(smp_s2s)];
      sts(end+1)=min(s1,s2s);
    else
      sts(end+1) = min(zj);
    end
   
    idx_e1=setdiff(1:size(zv,2),smp_e1);
    if ~isempty(idx_e1)
      zv_e1=zv(:,idx_e1);
      [s2e,e2e,smp_s2e,smp_e2e]=add_wide_peak(zv_e1,0,1);
      sample_e=[smp_e1 idx_e1(smp_e2e)];
      ens(end+1)=max(e1,e2e);
    else
      ens(end+1) = max(zj);
    end
  else
    error('k not supproted');
  end
else
  if k==-1
    for i=1:(length(samples))
      %    zsc=(zsc_c+sum(zv(:,setdiff(1:size(zv,2),i)),2));
      zsc=zsctot-zv(:,i);
      [mx,mi]=max(zsc);
      rg=find(zsc==mx);
      sts(end+1)=min(rg);
      ens(end+1)=max(rg);
    end
  elseif k==-2
    for i=1:(length(samples)-1)
      for j=(i+1):length(samples)
        %      zsc=(zsc_c+sum(zv(:,setdiff(1:size(zv,2),[i j])),2));
        zsc=zsctot-sum(zv(:,[ i j]),2);
        [mx,mi]=max(zsc);
        rg=find(zsc==mx);
        sts(end+1)=min(rg);
        ens(end+1)=max(rg);
      end
    end
  elseif k==-3
    for i=1:(length(samples)-2)
      for j=(i+1):(length(samples)-1)
        for n=(j+1):length(samples)
          %        zsc=(zsc_c+sum(zv(:,setdiff(1:size(zv,2),[i j n])),2));
          zsc=zsctot-sum(zv(:,[ i j n]),2);
          [mx,mi]=max(zsc);
          rg=find(zsc==mx);
          sts(end+1)=min(rg);
          ens(end+1)=max(rg);
        end
      end
    end
  else
    error('k>3 not supported');
  end
end

s=min(find(zj==min(sts)))+chr_zero;
e=max(find(zj==max(ens)))+chr_zero;
