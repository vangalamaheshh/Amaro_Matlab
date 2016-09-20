function Ps = enumerate_partitions(n,k)
%  Efficiently generates all numerical partitions of n 
  % k is an optional parameter.  If k is present, the program enumerates all
  % partitions with largest part less than or equal to k.
  % returns Ps, a cell array containing the partitions.
  
  % This is a matlab translation of C program downloaded from Frank Ruskey's
  % Combinatorial Object Server
    
    
    global p Ps idx
    Ps = {};
    idx = 1;
    if ~exist('k','var') || isempty(k) || k >= n
      % in this case, find all partitions
      P(2*n,n,0);
    else
      for j=1:k
        P(n,j,1);
      end
    end

function P(n,k,t)
  global p Ps idx
  p(t+1) = k;
  if n == k
    Ps{idx} = p(2:t+1)';
    idx = idx+1;
  end
  for j=min(k,n-k):-1:1
    P(n-k,j,t+1);
  end
    
    
    
      
