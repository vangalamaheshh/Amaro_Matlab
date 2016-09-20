function [D,use]=remove_D_types(D,typec,use_names)

dontuse=[];
for i=1:length(typec)
  if use_names
    dontuse=union(dontuse,findstrings(D.sdesc,typec{i}));
  else
    pos=findstrings(D.supdesc,typec{i});
    if length(pos)==1
      dontuse=union(dontuse,find(D.supdat(pos,:)));
    else
      disp(['Non unique identifier:' typec{i}]);
    end
  end
end
use=setdiff(1:size(D.dat,2),dontuse);

D=reorder_D_cols(D,use);

