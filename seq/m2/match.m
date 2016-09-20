function [ p_value ] = match( hist, norm_hist, mut_obs )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

metric = min(diff(sort(mut_obs)));
%keyboard
p_value = sum(norm_hist(1:metric));
    

end

