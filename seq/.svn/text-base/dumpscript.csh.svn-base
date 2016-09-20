#! /bin/csh -f

if ($#argv != 3 & $#argv != 4) then
  echo "Usage: dumpscript bamfile baifile outfile [banner]"
  exit 1
endif

set bamfile = $1
set baifile = $2
set outfile = $3

echo "Step 1: DumpReads (java)"
java -classpath /xchip/tcga/gbm/analysis/lawrence/samtools/trunk/sam-jdk/dist/sam-1.0.jar:/xchip/tcga/gbm/analysis/lawrence/sam DumpReads $bamfile $baifile $outfile.1

echo "Step 2: sort (unix)"
mkdir $outfile.tmp
sort $outfile.1 -T $outfile.tmp > $outfile.2
rmdir $outfile.tmp

echo "Step 3: joindump (perl)"
cat $outfile.2 | perl /xchip/tcga/gbm/analysis/lawrence/sam/joindump.pl > $outfile.3

mv $outfile.3 $outfile
rm $outfile.1 $outfile.2
echo "Done!"

