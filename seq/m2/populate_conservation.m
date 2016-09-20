function [ conservation_track ] = populate_conservation( genomic_array, fileCons, chr )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if isempty(genomic_array) || isempty(chr)
  fprintf('populate_conservation: WARNING empty input, empty output\n');
  conservation_track = [];
else
  try
     chr = repmat(chr, length(genomic_array), 1);
     conservation_track = fileCons.get(chr, genomic_array);
     conservation_track = (double(conservation_track)-50)/(25/3);      % convert back to negative to positive
     conservation_track(conservation_track==18) = NaN; % missing data
  catch me
    disp(me);
    disp(me.message);
    fprintf('ERROR IN populate_conservation\n');
    keyboard
  end
end
