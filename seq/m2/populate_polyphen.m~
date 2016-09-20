function [ polyphen_track ] = populate_polyphen( genomic_array, fA, fC, fG, fT, chr)
%Constructs polyphen track. 

    polyphen_track = zeros(size(genomic_array, 2), 4); 
    for i = 1:size(genomic_array, 2)
        polyphen_track(i, 1) = fA.get(chr, genomic_array(i));
        polyphen_track(i, 2) = fC.get(chr, genomic_array(i));
        polyphen_track(i, 3) = fG.get(chr, genomic_array(i));
        polyphen_track(i, 4) = fT.get(chr, genomic_array(i));
  %      keyboard
    end
  %  disp(polyphen_track)
%    keyboardx
end

