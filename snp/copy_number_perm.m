function [h_amp,h_del,naamp,nadel]=copy_number_perm(Csmooth,niter,method)

if ischar(method)
  method.method=method;
end

if niter==0
  noperm=1;
  niter=1;
else
  noperm=0;
end

if isfield(method,'t1s')
  t1s=method.t1s;
  t2s=method.t2s;
else
  t1s=NaN;
  t2s=NaN;
end

h_amp=zeros(1001,length(t1s));
h_del=zeros(1001,length(t2s));

for j=1:niter
  dat=Csmooth;
  if ~noperm
    for i=1:size(dat,2)
      r=randperm(size(dat,1));
      dat(:,i)=dat(r,i);
    end
  end

  for k=1:length(t1s)
    t1=t1s(k);
    t2=t2s(k);
    switch method.method
     case 'min2'
      if isnan(t1)
        t1=4/3; % 1.6;
        t2=3; % 2.5;
      end
      %  del=dat<t1;
      amp=dat>t2;
      
      namp=sum(amp,2);
      %  ndel=sum(del,2);
      %xamp=Corig.dat;
      xamp=dat;
      xamp(~amp)=NaN;
      aamp=nanmean(xamp,2);
      naamp=namp.*(aamp-2);
     case 'log'
      % t1=log2(1.6)-1; % 4/3, 1.6;
      % t2=log2(2.5)-1; % 3, 2.5;
      
      % t1=log2(4/3)-1; % 4/3, 1.6;
      % t2=log2(3)-1; % 3, 2.5;
      
      % t1=-0.4;
      % t2=0.4;
      
      if isnan(t1)
        cbs=1;
        if cbs
          %      t1=log2(1.68)-1;
          %      t2=log2(2.4)-1;
          t1=-0.25;
          t2=0.25;
        else
          t1=-0.6;
          t2=0.6;
        end 
      end

      del=dat<t1;
      amp=dat>t2;
      namp=sum(amp,2);
      ndel=sum(del,2);
      %xamp=Corig.dat;
      xamp=dat;
      xamp(~amp)=NaN;
      
      aamp=nanmean(xamp,2);
      naamp=namp.*aamp;
      
      xdel=dat;
      xdel(~del)=NaN;
      adel=nanmean(xdel,2);
      nadel=-ndel.*adel;
     
     case 'n'
      if isnan(t1)
        cbs=1;
        if cbs
          t1=-0.25;
          t2=0.25;
        else
          t1=-0.6;
          t2=0.6;
        end 
      end

      del=dat<t1;
      amp=dat>t2;
      namp=sum(amp,2);
      ndel=sum(del,2);
      
      naamp=namp;
      nadel=ndel;

    end
    %xdel=Corig.dat;
    %  xdel=C.dat;
    %  xdel(~C.del)=NaN;
    %  C.adel=nanmean(xdel,2);
    %  y=round(naamp*10)/10;
    y=naamp;
    y(isnan(y))=0;
    h_amp(:,k)=h_amp(:,k)+histc(y,0:0.1:100);
    
    %  y=round(nadel*10)/10;
    y=nadel;
    y(isnan(y))=0;
    h_del(:,k)=h_del(:,k)+histc(y,0:0.1:100);
  end

  
  disp(j);
end


if (0)
  C.dat=C.smooth;
  t1=4/3; % 1.6;
  t2=3; % 2.5;
        %  C.del=C.dat<t1;
  C.amp=C.dat>t2;
  
  C.namp=sum(C.amp,2);
  %  C.ndel=sum(C.del,2);
  %xamp=Corig.dat;
  xamp=C.dat;
  xamp(~C.amp)=NaN;
  C.aamp=nanmean(xamp,2);
  %xdel=Corig.dat;
  %  xdel=C.dat;
  %  xdel(~C.del)=NaN;
  %  C.adel=nanmean(xdel,2);
  C.naamp=C.namp.*(C.aamp-2);
  y=round(C.naamp*10)/10;
  y(isnan(y))=0;
end
