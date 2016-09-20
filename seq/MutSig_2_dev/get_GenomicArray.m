function [ genomic_array ] = get_GenomicArray( exons, length )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
genomic_array = zeros(1, length); 
counter = 1;

for i = 1:size(exons,1)
     current = exons(i,1);
     while current <= exons(i, 2)
         genomic_array(counter) = current;
         current = current + 1;
         counter = counter + 1;
 %        keyboard
     end 
end

