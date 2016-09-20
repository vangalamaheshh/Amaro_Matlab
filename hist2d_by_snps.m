function h_counts = hist2d_by_snps(M,length_vector,x_edges,y_edges)  
% ---
% $Id$
% $Date: 2008-05-15 10:45:50 -0400 (Thu, 15 May 2008) $
% $LastChangedBy: cmermel $
% $Rev$
  
  %% hist2d_by_snps counts the number of elements in the Nx2 matrix M =
  %%(x,y) residing in each bin defined by the two edge vectors
  %%x_edges (a p x 1 vector) and y_edges (a q x 1 vector)  where each
  %%element has a multiplicity given by the N x1 length_vector
    
  %% The function returns a pxq matrix h_counts, where h_counts(i,j) is
  %%the sum of the number of elements with x-value between x_edges(i) and
  %%x_edges(i+1) and y-value between y_edges(j) and y_edges(j+1).
    
  h_counts = zeros(length(x_edges),length(y_edges));
  
  x = M(:,1);
  y = M(:,2);
  
  for i=1:length(x_edges)-1
    cur_x = intersect(find(x>x_edges(i)),find(x<=x_edges(i+1)));
    cur_y = y(cur_x);
    cur_lengths = length_vector(cur_x);
    for j=1:length(y_edges)-1
      h_counts(i,j) = sum(cur_lengths(intersect(find(cur_y>y_edges(j)), ...
                                                find(cur_y<=y_edges(j+1)))));
    end
    
    h_counts(i,j+1) = sum(cur_lengths(find(cur_y>y_edges(end))));
  end
  
  %% Fill in last-column (i=length(x_edges))
  
  cur_x = find(x>x_edges(end));
  cur_y = y(cur_x);
  cur_lengths = length_vector(cur_x);
  
  for j=1:length(y_edges)-1
    h_counts(length(x_edges),j) = sum(cur_lengths(intersect(find(cur_y> ...
                                                      y_edges(j)),find(cur_y<=y_edges(j+1)))));
  end
  h_counts(length(x_edges),j+1) = sum(cur_lengths(find(cur_y> ...
                                                    y_edges(end))));
  
  
  
  
      
      
    
    
    