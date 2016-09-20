function [ p, index, ci_ratio, ci_bounds, k, n ] = p_value(h, observed,  bins)

    [tmpzz,index]=min(abs(bins-observed));
    n = sum(h);
    k = sum(h(index:size(bins,2)));

    [p, ci_bounds] = binofit(k, n);
    
    ci_ratio = 2*1.96*sqrt((n-k+1) / ((k+1)* (n+3)));

end

