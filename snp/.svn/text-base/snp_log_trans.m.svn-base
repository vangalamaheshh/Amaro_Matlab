function D=snp_log_trans(D,rev,force)

thresh=0.1;
perf=0;
if exist('rev','var') && rev==1
  if abs(mean(D.dat(:)))>0.5
    warning('Seems like non-log data.');
    if  exist('force','var') && force==1
      disp('Performing anyway');
      perf=1;
    else
      disp('Ignoring');
    end
  else
    perf=1;
  end
  if perf
    D.dat=2.^(D.dat+1);
  end
else
  if abs(mean(D.dat(:)))<0.5 
    warning('Seems like already log-transformed data. Check!');
    if  exist('force','var') && force==1
      disp('Performing anyway');
      perf=1;
    else
      disp('Ignoring');
    end
  else
    perf=1;
  end
  if perf
    D.dat(D.dat<thresh)=thresh;
    D.dat=log2(D.dat)-1;
  end
end

