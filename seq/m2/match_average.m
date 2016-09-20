function [ p_mle] = match_average(ks_stat_array, mutations, length, mean_curve )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    mut_curve = calculate_curve(mutations, length);
%    keyboard
    mut_ks_stat = max(as_row(mut_curve) - mean_curve');
    n = nnz(ks_stat_array); 
    k = size(find(ks_stat_array >= mut_ks_stat), 1);
    p_mle = k/n;
 %   keyboard


end

