function zmean = mean_profile_over_arms(z,C,cyto,fract_n)
% For a score z spanning the genome, this function determines
% the average profile of z across a chromosome arm. For chromosome
% arms with greater than fract_n markers, the markers are equally
% distributed among fract_n bins from telomere to centromere, and the
% mean across all these markers is determined.
% Chromosome arms are weighted by the number of markers they contain.
%
% ---
% $Id$
% $Date: 2009-08-14 10:00:01 -0400 (Fri, 14 Aug 2009) $
% $LastChangedBy: rameen $
% $Rev$

if isfield(C,'armn')
  D=C;
  clear C;
else
  D=add_cyto(C,cyto);
end

if size(z,1)==1
  z = z';
end

zsum = zeros(fract_n,1);
counter = zsum;
zk = zsum;
zz=z;

for i = 1:max(D.chrn)
  for j = 1:2
    inchrarm = find(D.chrn==i & D.armn ==j);
    if length(inchrarm)>fract_n
      if j==2
        zz(inchrarm)=flipud(z(inchrarm));
       end
      iarm = ceil([1:length(inchrarm)]*fract_n/length(inchrarm));
      for k = 1:fract_n
        zmarkers=find(iarm==k);
        zk(k) = sum(zz(inchrarm(zmarkers)));
        counter(k) = counter(k) + length(zmarkers);
      end
      zsum=zsum+zk;
    end
  end
end

if min(counter) > 0
  zmean = zsum./counter;
else
  error('No chromosome arms had more than %f markers',fract_n)
end
