function h_counts = hist_by_snps(y,num_snp_vector,edges)
% ---
% $Id$
% $Date: 2007-11-13 11:39:52 -0500 (Tue, 13 Nov 2007) $
% $LastChangedBy: cmermel $
% $Rev$
  
  for i=1:(length(edges)-1)
    h_counts(i) = sum(num_snp_vector(intersect(find(y>edges(i)),find(y<edges(i+1)))));
  end
  
  h_counts(i+1) = sum(num_snp_vector(find(y>edges(end))));
  
 
  
  
