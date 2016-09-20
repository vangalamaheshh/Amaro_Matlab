function jobid = fishall(sample,span,banner)

R = dRanger_load_results(sample);
R = reorder_struct(R,R.normreads==0 & R.filterHCLB==0);
nr = slength(R);
fprintf('fishall:  %s (%d somatic RAs) --> "%s"\n',sample,nr,banner);
R.bait = cell(nr,1);
s1 = ceil(span/2); s2 = span-s1;
for i=1:nr
  span1 = R.max1(i)-R.min1(i)+1;
  span2 = R.max2(i)-R.min2(i)+1;
  fprintf('%d %d %d (%d)\n',i,span1,span2,R.tumreads(i));
  if R.str1(i)==0, part1 = genome_region(R.chr1(i),R.max1(i)-s1,R.max1(i)+s2);
  else part1 = rc(genome_region(R.chr1(i),R.min1(i)-s1,R.min1(i)+s2)); end
  if R.str2(i)==0, part2 = genome_region(R.chr2(i),R.max2(i)-s1,R.max2(i)+s2);
  else part2 = rc(genome_region(R.chr2(i),R.min2(i)-s1,R.min2(i)+s2)); end
  R.bait{i} = [num2str(i) ' ' part1 ' 12x12 -12:76 ' part2 ' 12x12'];
end
fishdir = ['/xchip/tcga_scratch/lawrence/' sample '/fishall'];
mkdir(fishdir); save_struct(R,[fishdir '/somatic_RA_list.txt']);
jobid = fishreads(sample,'tumor',R.bait,'fishall',banner);
