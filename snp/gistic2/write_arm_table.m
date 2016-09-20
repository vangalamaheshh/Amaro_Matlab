function write_arm_table(fname,A,cyto)
%WRITE_ARM_TABLE save arm-level copy number data to file
%
%   WRITE_ARM_TABLE(FNAME,A,CYTO,RG)
%
% Outputs A (copy number reduced to arm level) to the file FNAME, 
% using CYTO to get arm names.

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  nsamples = size(A.dat,2);
  n_arms = size(A.dat,1);
  fid = fopen(fname,'w');

  % print header
  fprintf(fid,'Arm');
  for j=1:nsamples
    fprintf(fid,'\t%s',char(A.sdesc(j,:)));
  end
  fprintf(fid,'\n');
  
  fmt = ['%s' repmat('\t%1.3f',1,nsamples) '\n'];
  for j = 1:n_arms
      band = regexp(cyto(A.cyton(j)).name,'([0-9X]+[pq])','tokens');
      fprintf(fid,fmt,band{1}{1},A.dat(j,:)); 
  end  
  fclose(fid);

