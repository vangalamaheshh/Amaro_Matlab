function [ categories_track ] = populate_categories( genomic_array, fileCat, chr )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    categories_track = zeros(size(genomic_array, 2), 1); 
    
    for i = 1:size(genomic_array, 2)
        categories_track(i, 1) = fileCat.get(chr, genomic_array(i));

  %      keyboard
    end

end

