function s=display_soap_result(x,nt)

if ~exist('nt','var')
  nt=0;
end

s=[];
tab=char(9);
nl=sprintf(newline);
if isempty(x)
  s=[s '---'];
elseif isstruct(x)
  s=[s sprintf(['{' newline])];
  if nt>0
    s=[s repmat(tab,1,nt)];
  end
  for i=1:length(x)
    f=fieldnames(x(i));
    for j=1:length(f)
      s=[s sprintf([f{j} ':'])];
      s=[s tab];
      s=[s display_soap_result(getfield(x,{i},f{j}),nt+1) nl];
      if nt>0
        s=[s repmat(tab,1,nt)];
      end
    end
  end
  s=[s sprintf(['}'])];
elseif iscell(x)
  for i=1:length(x)
    if i<length(x)
      s=[s display_soap_result(x{i},nt+1) nl];
      if nt>0
        s=[s repmat(tab,1,nt)];
      end
    else
      s=[s display_soap_result(x{i},nt+1)];
    end
  end
elseif ischar(x)
  s=[s x];
elseif isnumeric
  s=[s sprintf('%f',x)];
end


if nargout==0
  fprintf(1,'%s',s);
end
