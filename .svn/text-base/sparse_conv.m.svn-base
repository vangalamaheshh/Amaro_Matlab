function x = sparse_conv(a,b,max_sparsity)
  
  %% This function implements the convolution of two (column) vectors a
  %% and b, a*b, supporting the possibility that one or both of a and b
  %% are represented as sparse matrices.
    
    if ~exist('max_sparsity','var') || isempty('max_sparsity')
      max_sparsity = .02;
    end
    
    if size(a,1) == 1 && size(a,2) > 1
      a = a';
    end
    
    if size(b,1) == 1 && size(b,2) > 1
      b = b';
    end
    
    if size(a,2) > 1 || size(b,2) > 1
      error('inputs must be a row or column vector')
    end
    
    if size(a,1) > size(b,1) 
      %% Want a to be shorter vector, so swap a and b
      temp_a = a;
      a = b;
      b = temp_a;
    end
      
    if ~issparse(a) && ~issparse(b)
      x = conv(a,b);
    elseif max(length(find(a))/length(a),length(find(b))/length(b)) >= ...
          max_sparsity
      x = conv(full(a),full(b));
    else
      x = sparse(length(a)+length(b)-1,1);
      [ai,aj,ax] = find(a);
      [bi,bj,bx] = find(b);
      for i=1:length(ai)
        for j=1:length(bi)
          x(ai(i)+bi(j)-1) = x(ai(i)+bi(j)-1)+a(ai(i))*b(bi(j));
        end
      end
    end
    
    
      
