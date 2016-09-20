function [ curve] = calculate_curve_fast( mutations, gene_length )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%matrix = zeros(length(mutations)*(length(mutations)-1), 1);
%counter = 1; 
bins = 0:gene_length - 1;
histogram = zeros(length(bins),1);
%pre_array = zeros((length(mutations)*(length(mutations)-1))/2,1);
for i=1:length(mutations)-1
     for j = i+1:length(mutations)
 %        pre_array(counter) = abs(mutations(i) - mutations(j)); 
          difference = abs(mutations(i) - mutations(j));
          histogram(difference+1) = histogram(difference+1) + 1;
 %        counter = counter + 1;
     end
end
%disp(pre_array);
%keyboard
curve = cumsum(histogram);
%keyboard




