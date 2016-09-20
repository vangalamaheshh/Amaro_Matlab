function [ coverage_track ] = populate_coverage( genomic_array, fileCov, chr )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    coverage_track = zeros(size(genomic_array, 2), 1); 
    %keyboard
    for i = 1:size(genomic_array, 2)
        coverage_track(i, 1) = fileCov.get(chr, genomic_array(i));

  %      keyboard
    end
 %   disp(polyphen_track)
end

