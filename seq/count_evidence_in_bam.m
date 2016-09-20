function [ref mut other1 other2] = count_evidence_in_bam(M,bam,first_idx,last_idx,blacklist)
% Mike Lawrence 2009-11-05

if ~exist('blacklist','var'), blacklist = 'none'; end

require_fields(M,{'chr','start','ref_allele','tum_allele1','tum_allele2'});

M.chrn = convert_chr(M.chr);
M = make_numeric(M,{'start'});

if ~exist('first_idx','var'),first_idx=1;end
if ~exist('last_idx','var'),last_idx=slength(M);end

base('ACGTN')=1:5;

javaclasspath('/xchip/tcga/gbm/analysis/lawrence/samtools/sam-1.07.jar',...
   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq');
import net.sf.samtools.*;
import java.io.*;
import java.lang.*;

demand_file(bam);
x = BamGrasp();
x.setQuietModeOn;
x.openFile(String(bam),String(blacklist));

try

ref = nan(slength(M),1);
mut = nan(slength(M),1);
other1 = nan(slength(M),1);
other2 = nan(slength(M),1);
for i=first_idx:last_idx, fprintf('%d/%d ',i,slength(M));
  [R B S] = BamGrasp_load_region(x,M.chrn(i),M.start(i));
  if strcmp(M.ref_allele{i},M.tum_allele1{i}), newbase = M.tum_allele2{i}; else newbase = M.tum_allele1{i}; end
  bidx = find(B(:,4)==M.start(i));
  ref(i) = sum(B(bidx,1)==base(S));
  mut(i) = sum(B(bidx,1)==base(newbase));
  others = setdiff(1:4,base([newbase S]));
  other1(i) = sum(B(bidx,1)==others(1));
  other2(i) = sum(B(bidx,1)==others(2));
end

x.closeFile();

catch me,excuse(me);end




