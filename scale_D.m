function [D,scale_value]=scale_D(D,scaling)

D=add_history(D,mfilename,scaling);

if ischar(scaling)
  scaling1.method=scaling;
  scaling=scaling1;
end

if ~isfield(scaling,'value')
  scaling.value=1;
  scaling.use_median=1;
end

switch scaling.method
 case 'mean'
  sf=repmat(scaling.value,1,size(D.dat,2))./nanmean(D.dat,1);
  if isfield(scaling,'use_median') && scaling.use_median
    [sf,scale_value]=use_mean(sf,scaling.value);
  end
  D.dat=D.dat.*repmat(sf,size(D.dat,1),1);
  D.sscale=sf;
 case 'median'
  sf=repmat(scaling.value,1,size(D.dat,2))./nanmedian(D.dat,1);
  if isfield(scaling,'use_median') && scaling.use_median
    [sf,scale_value]=use_median(sf,scaling.value);
  end
  D.dat=D.dat.*repmat(sf,size(D.dat,1),1);
  D.sscale=sf; 
 case 'prctile'
  sf=repmat(scaling.value,1,size(D.dat,2))./prctile(D.dat,scaling.prctile);
  if isfield(scaling,'use_median') && scaling.use_median
    [sf,scale_value]=use_median(sf,scaling.value);
  end
  D.dat=D.dat.*repmat(sf,size(D.dat,1),1);
  D.sscale=sf;
 case 'meanoftop'
  X=sort(D.dat,1);
  sf=repmat(scaling.value,1,size(D.dat,2))./nanmean(X(round((1-scaling.fraction)*size(X,1)):end,:),1);
  if isfield(scaling,'use_median') && scaling.use_median
    [sf,scale_value]=use_median(sf,scaling.value);
  end
  D.dat=D.dat.*repmat(sf,size(D.dat,1),1);
  D.sscale=sf;
 otherwise
  error('no such scaling method');
end

D.sscale=num2str(D.sscale');

function [sf,scaling_value]=use_median(sf,sv)
medsf=median(sf);
scaling_value1=sv/medsf;
scaling_value=roundsig(scaling_value1,2);
sf=sf/(sv/scaling_value);
