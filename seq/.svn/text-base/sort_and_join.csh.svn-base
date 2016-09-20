#! /bin/csh -f

if ($#argv != 2 & $#argv != 3) then
  echo "Usage: sort_and_join infile outfile [banner]"
  exit 1
endif

set infile = $1
set outfile = $2

echo "Step 1: sort (unix)"
mkdir $outfile.tmp
sort $infile -T $outfile.tmp > $outfile.1
rmdir $outfile.tmp

echo "Step 2: joindump (perl)"
cat $outfile.1 | perl /xchip/tcga/gbm/analysis/lawrence/sam/joindump.pl > $outfile.2

mv $outfile.2 $outfile
rm $outfile.1
echo "Done!"

