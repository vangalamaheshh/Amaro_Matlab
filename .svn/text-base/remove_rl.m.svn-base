function rl=remove_rl(rl,rli,v)

if iscell(rl)
  for i=1:length(rl)
    if ~isempty(rli{i})
      error('do not use!!');
      % take care of changing indices
      rl{i}=remove_rl(rl{i},rli{i});
    end
  end
else
  if rli==1
    rl=rl(2:end,:);
  elseif rli==size(rl,1)
    rl=rl(1:(end-1),:);
  else
    if (rl(rli-1,3)==rl(rli+1,3)) && ...
      ((exist('v','var') && v(rl(rli-1,1))==v(rl(rli+1,1))) || (~exist('v','var')))
      en=rl(rli+1,2);
      rl=rl(setdiff(1:size(rl,1),[rli rli+1]),:);
      rl(rli-1,2)=en;
    else
      rl=rl(setdiff(1:size(rl,1),rli),:);      
    end
  end
end
