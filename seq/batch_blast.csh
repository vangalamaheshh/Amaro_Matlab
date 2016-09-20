#! /bin/csh -f

if ($#argv != 2) then
  echo "Usage: batch_blast <qfile> <hfile>"
  exit 1
endif

set qfile = $1
set hfile = $2

set dbfile = "/xchip/tcga/gbm/analysis/lawrence/genome/hg18/orig/hg18.fa"

blastall -p blastn -d $dbfile -m8 -e1e-50 -nT -FF < $qfile > $hfile.partial

if ($status == 0) then
  mv $hfile.partial $hfile
endif
