function [new_H new_h] = calc_range_dist(d,n)
  
  % Function takes a distribution of background g-scores (d) and computes
  % the range distribution h for a sample of size n from that
  % distribution. In other words, it calculates the distribution of r(n), the
  % difference between the sample max and sample min in a sample of size
  % n.
    
% H = cumulative distribution of range
  % h = distribution of range  
    
  if size(d,2) == 1
    new_d{1} = d;
    h = {zeros(size(d,1),1)};
    H = h;
  else
    new_d = d;
    h = {zeros(size(d{1},1),1), zeros(size(d{1},1),1)};
    H = h;
  end 
   
    if ~exist('d','var') | isempty(d)
      error('Must supply snp_score background distribution!')
    end
    
    if ~exist('n','var') | isempty(n)
      n=1;
    end
    
    for k=1:size(new_d,2)
      D{k} = cumsum(new_d{k});
      for r=0:size(new_d{k},1)-1
        %if mod(r,1000) == 0
        %  disp(['r = ' num2str(r)]);
        %end
        h{k}(r+1) = 0;
        for x=r+1:size(new_d{k},1)
          h{k}(r+1) = h{k}(r+1) + n*(n-1)*new_d{k}(x)*new_d{k}(x-r)*(D{k}(x)-D{k}(x- ...
                                                            r))^(n-2);
        end
      end
      h{k} = h{k}/sum(h{k});
      H{k} = cumsum(h{k});
    end
     
    if size(h,1) == 1
      new_h = h{1};
      new_H = H{1};
    else
      new_h = h;
      new_H = H;
    end
    
    
