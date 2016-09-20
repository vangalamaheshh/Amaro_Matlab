function write_arm_medians(fname,D,arm_medians,names)
  

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  verbose('Writing arm medians file...',30);
  
  f = fopen(fname,'w');
  fprintf(f,'%s','Chromosome Arm');
  for j=1:size(arm_medians,2)
    fprintf(f,'\t%s',char(D.sdesc(j,:)));
  end
  fprintf(f,'\n');
    
  for j=1:length(names)
    fprintf(f,'%s',char(names(j)));
    for i=1:size(arm_medians,2)
      fprintf(f,'\t%1.3f',arm_medians(j,i));
    end
    fprintf(f,'\n');
  end

  
