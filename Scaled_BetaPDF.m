function PDF = Scaled_BetaPDF(y, a, b, p, q)
PDF = ( (y-p).^(a-1) .* (q - y).^(b-1) ) ./ ( (q - p).^(a+b-1) .* beta(a,b) );