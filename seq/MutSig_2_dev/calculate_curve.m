function [ curve] = calculate_curve( mutations, gene_length )
%%% OBSOLETE
%%% replaced by fast_calculate_curve.c

%method=1;

%if method==1

  histogram = zeros(gene_length,1);
  for i=1:length(mutations)-1
    for j = i+1:length(mutations)
      difference = abs(mutations(i) - mutations(j));
      histogram(difference+1) = histogram(difference+1) + 1; 
    end 
  end

%elseif method==2       % ~60% SLOWER!

%  difference = abs(bsxfun(@minus,mutations,mutations'));
%  histogram = histc(difference(:),0:(gene_length-1));
%  histogram(1) = histogram(1) - length(mutations); % subtract diagonal

%end

curve = cumsum(histogram);
curve = curve / curve(end);
