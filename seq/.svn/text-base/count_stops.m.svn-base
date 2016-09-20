function stops = count_stops(dna,frame)

% frame 1,2,3 = strand as given, +0, +1, +2
% frame 4,5,6 = revcomp strand, +0, +1, +2

if frame==1
   i=1;
elseif frame==2
   i=2;
elseif frame==3
   i=3;
elseif frame==4
   dna = my_seqrcomplement(dna);
   i=1;
elseif frame==5
   dna = my_seqrcomplement(dna);
   i=2;
elseif frame==6
   dna = my_seqrcomplement(dna);
   i=3;
else
   error('invalid frame %d\n', frame);
end

n = length(dna);
dna = upper(dna);

stops = 0;
while(i+2 <= n)
   codon = dna(i:i+2);
   if strcmp(codon, 'TAA') || strcmp(codon, 'TAG') ...
      || strcmp(codon, 'TGA')
         stops = stops + 1;
   end
   i=i+3;
end

end
